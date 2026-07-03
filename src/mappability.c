#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "mbpriv.h"
#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, 0x10000)

typedef struct {
	int32_t st, en, max_en;
	uint8_t cap;
} mb_mapq_intv_t;

typedef struct {
	int32_t n, m;
	mb_mapq_intv_t *a;
} mb_mapq_intv_v;

struct mb_mapq_track_s {
	int32_t n_ctg;
	mb_mapq_intv_v *ctg;
};

static int cmp_intv(const void *a, const void *b)
{
	const mb_mapq_intv_t *x = (const mb_mapq_intv_t*)a;
	const mb_mapq_intv_t *y = (const mb_mapq_intv_t*)b;
	if (x->st != y->st) return x->st < y->st? -1 : 1;
	if (x->en != y->en) return x->en < y->en? -1 : 1;
	return (int)x->cap - (int)y->cap;
}

static int32_t clamp_cap(int32_t cap)
{
	if (cap < 0) return 0;
	if (cap > 60) return 60;
	return cap;
}

static int32_t cap_from_score(const char *s, int32_t default_cap)
{
	char *end = 0;
	double v;
	if (s == 0 || *s == 0) return clamp_cap(default_cap);
	errno = 0;
	v = strtod(s, &end);
	if (s == end || errno == ERANGE || isnan(v)) return clamp_cap(default_cap);
	if (v <= 0.0) return 0;
	if (v <= 1.0) return clamp_cap((int32_t)(60.0 * v + .499));
	if (v <= 60.0) return clamp_cap((int32_t)(v + .499));
	if (v <= 100.0) return clamp_cap((int32_t)(60.0 * v / 100.0 + .499));
	return 60;
}

static char *next_field(char **p_)
{
	char *p = *p_, *q;
	while (*p && isspace((unsigned char)*p)) ++p;
	if (*p == 0) {
		*p_ = p;
		return 0;
	}
	q = p;
	while (*q && !isspace((unsigned char)*q)) ++q;
	if (*q) *q++ = 0;
	*p_ = q;
	return p;
}

static int name_eq_alias(const char *a, const char *b)
{
	if (strcmp(a, b) == 0) return 1;
	if (strncmp(a, "chr", 3) == 0 && strcmp(a + 3, b) == 0) return 1;
	if (strncmp(b, "chr", 3) == 0 && strcmp(a, b + 3) == 0) return 1;
	if ((strcmp(a, "M") == 0 && strcmp(b, "MT") == 0) || (strcmp(a, "MT") == 0 && strcmp(b, "M") == 0)) return 1;
	if (strncmp(a, "chr", 3) == 0 && strncmp(b, "chr", 3) == 0) {
		const char *x = a + 3, *y = b + 3;
		return (strcmp(x, "M") == 0 && strcmp(y, "MT") == 0) || (strcmp(x, "MT") == 0 && strcmp(y, "M") == 0);
	}
	return 0;
}

static int32_t find_tid(const l2b_t *l2b, const char *name, int32_t last_tid)
{
	int32_t i;
	if (last_tid >= 0 && last_tid < (int32_t)l2b->n_ctg && name_eq_alias(l2b->ctg[last_tid].name, name))
		return last_tid;
	for (i = 0; i < (int32_t)l2b->n_ctg; ++i)
		if (name_eq_alias(l2b->ctg[i].name, name))
			return i;
	return -1;
}

static void add_intv(mb_mapq_intv_v *v, int32_t st, int32_t en, int32_t cap)
{
	mb_mapq_intv_t *last;
	if (cap >= 60 || st >= en) return;
	if (v->n > 0) {
		last = &v->a[v->n - 1];
		if (last->en >= st && last->cap == cap) {
			if (last->en < en) last->en = en;
			return;
		}
	}
	kom_grow(mb_mapq_intv_t, v->a, v->n, v->m);
	v->a[v->n].st = st;
	v->a[v->n].en = en;
	v->a[v->n].cap = (uint8_t)cap;
	++v->n;
}

static void sort_and_merge(mb_mapq_intv_v *v)
{
	int32_t i, k;
	if (v->n == 0) return;
	if (v->n == 1) {
		v->a[0].max_en = v->a[0].en;
		return;
	}
	qsort(v->a, v->n, sizeof(*v->a), cmp_intv);
	for (i = k = 0; i < v->n; ++i) {
		if (k > 0 && v->a[k - 1].en >= v->a[i].st && v->a[k - 1].cap == v->a[i].cap) {
			if (v->a[k - 1].en < v->a[i].en) v->a[k - 1].en = v->a[i].en;
		} else {
			if (k < i) v->a[k] = v->a[i];
			++k;
		}
	}
	v->n = k;
	for (i = 0; i < v->n; ++i)
		v->a[i].max_en = (i > 0 && v->a[i - 1].max_en > v->a[i].en)? v->a[i - 1].max_en : v->a[i].en;
}

mb_mapq_track_t *mb_mapq_track_load(const mb_idx_t *idx, const char *fn, int32_t default_cap)
{
	const l2b_t *l2b = idx? idx->l2b : 0;
	mb_mapq_track_t *track = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t line = {0,0,0};
	int dret, last_tid = -1;
	int64_t n_line = 0, n_loaded = 0;

	if (l2b == 0 || fn == 0) return 0;
	fp = gzopen(fn, "rb");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	track = kom_calloc(mb_mapq_track_t, 1);
	track->n_ctg = (int32_t)l2b->n_ctg;
	track->ctg = kom_calloc(mb_mapq_intv_v, track->n_ctg);
	while (ks_getuntil(ks, KS_SEP_LINE, &line, &dret) >= 0) {
		char *p, *chrom, *st_s, *en_s, *score_s;
		uint64_t st0, en0;
		int32_t tid, st, en, cap;
		++n_line;
		p = line.s;
		while (*p && isspace((unsigned char)*p)) ++p;
		if (*p == 0 || *p == '#' || strncmp(p, "track", 5) == 0 || strncmp(p, "browser", 7) == 0) continue;
		chrom = next_field(&p);
		st_s = next_field(&p);
		en_s = next_field(&p);
		score_s = next_field(&p);
		if (chrom == 0 || st_s == 0 || en_s == 0) continue;
		tid = find_tid(l2b, chrom, last_tid);
		if (tid < 0) continue;
		last_tid = tid;
		st0 = strtoull(st_s, 0, 10);
		en0 = strtoull(en_s, 0, 10);
		if (en0 <= st0) continue;
		if (st0 >= l2b->ctg[tid].len) continue;
		if (en0 > l2b->ctg[tid].len) en0 = l2b->ctg[tid].len;
		st = st0 > INT32_MAX? INT32_MAX : (int32_t)st0;
		en = en0 > INT32_MAX? INT32_MAX : (int32_t)en0;
		cap = cap_from_score(score_s, default_cap);
		add_intv(&track->ctg[tid], st, en, cap);
		if (cap < 60) ++n_loaded;
	}
	free(line.s);
	ks_destroy(ks);
	gzclose(fp);
	for (int32_t i = 0; i < track->n_ctg; ++i)
		sort_and_merge(&track->ctg[i]);
	if (kom_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] loaded %ld mappability cap intervals from %ld BED rows\n",
				__func__, kom_realtime(), kom_percent_cpu(), (long)n_loaded, (long)n_line);
	return track;
}

void mb_mapq_track_destroy(mb_mapq_track_t *track)
{
	int32_t i;
	if (track == 0) return;
	for (i = 0; i < track->n_ctg; ++i)
		free(track->ctg[i].a);
	free(track->ctg);
	free(track);
}

static int32_t mb_mapq_track_cap(const mb_mapq_track_t *track, int64_t tid, int64_t st, int64_t en)
{
	const mb_mapq_intv_v *v;
	int32_t lo, hi, cap = 60;
	if (track == 0 || tid < 0 || tid >= track->n_ctg || st >= en) return 60;
	v = &track->ctg[tid];
	if (v->n == 0) return 60;
	lo = 0, hi = v->n;
	while (lo < hi) {
		int32_t mid = lo + ((hi - lo) >> 1);
		if (v->a[mid].st < en) lo = mid + 1;
		else hi = mid;
	}
	while (--lo >= 0) {
		if (v->a[lo].max_en <= st) break;
		if (v->a[lo].en > st && v->a[lo].cap < cap) cap = v->a[lo].cap;
	}
	return cap;
}

void mb_mapq_track_apply(const mb_mapq_track_t *track, int32_t n_regs, mb_hit_t *regs)
{
	int32_t i;
	if (track == 0) return;
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		int32_t cap;
		if (r->parent != r->id) continue;
		cap = mb_mapq_track_cap(track, r->tid, r->ts, r->te);
		if (r->mapq > cap) r->mapq = cap;
	}
}
