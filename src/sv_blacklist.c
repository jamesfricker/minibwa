#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mbpriv.h"

static int cmp_region(const void *a, const void *b)
{
	const mb_region_t *x = (const mb_region_t*)a, *y = (const mb_region_t*)b;
	if (x->st < y->st) return -1;
	if (x->st > y->st) return 1;
	if (x->en < y->en) return -1;
	if (x->en > y->en) return 1;
	return 0;
}

void mb_idx_clear_sv_blacklist(mb_idx_t *idx)
{
	int64_t i;
	if (idx == 0) return;
	if (idx->sv_bl) {
		for (i = 0; i < idx->l2b->n_ctg; ++i)
			free(idx->sv_bl[i]);
	}
	free(idx->sv_bl);
	free(idx->sv_bl_n);
	free(idx->sv_bl_m);
	idx->sv_bl = 0;
	idx->sv_bl_n = idx->sv_bl_m = 0;
	idx->n_sv_blacklist = 0;
}

static int find_ctg(const l2b_t *l2b, const char *name)
{
	int64_t i;
	const char *q = strncmp(name, "chr", 3) == 0? name + 3 : 0;
	for (i = 0; i < l2b->n_ctg; ++i) {
		const char *p = l2b->ctg[i].name;
		if (strcmp(p, name) == 0) return (int)i;
		if (q && strcmp(p, q) == 0) return (int)i;
		if (!q && strncmp(p, "chr", 3) == 0 && strcmp(p + 3, name) == 0) return (int)i;
	}
	return -1;
}

static void append_region(mb_idx_t *idx, int tid, int64_t st, int64_t en)
{
	int32_t n = idx->sv_bl_n[tid];
	kom_grow(mb_region_t, idx->sv_bl[tid], n, idx->sv_bl_m[tid]);
	idx->sv_bl[tid][n].st = st;
	idx->sv_bl[tid][n].en = en;
	idx->sv_bl_n[tid] = n + 1;
}

static void merge_regions(mb_idx_t *idx)
{
	int64_t tid;
	idx->n_sv_blacklist = 0;
	for (tid = 0; tid < idx->l2b->n_ctg; ++tid) {
		int32_t i, k, n = idx->sv_bl_n[tid];
		mb_region_t *a = idx->sv_bl[tid];
		if (n == 0) continue;
		qsort(a, n, sizeof(*a), cmp_region);
		for (i = k = 0; i < n; ++i) {
			if (k == 0 || a[i].st > a[k - 1].en) {
				a[k++] = a[i];
			} else if (a[k - 1].en < a[i].en) {
				a[k - 1].en = a[i].en;
			}
		}
		idx->sv_bl_n[tid] = k;
		idx->n_sv_blacklist += k;
	}
}

int64_t mb_idx_load_sv_blacklist(mb_idx_t *idx, const char *fn)
{
	FILE *fp;
	char line[4096];
	int64_t loaded = 0, skipped_ctg = 0, skipped_bad = 0;

	if (idx == 0 || fn == 0) return -1;
	fp = fopen(fn, "r");
	if (fp == 0) {
		if (kom_verbose >= 1)
			fprintf(stderr, "[ERROR] failed to open HMF SV blacklist BED '%s': %s\n", fn, strerror(errno));
		return -1;
	}

	mb_idx_clear_sv_blacklist(idx);
	idx->sv_bl_n = kom_calloc(int32_t, idx->l2b->n_ctg);
	idx->sv_bl_m = kom_calloc(int32_t, idx->l2b->n_ctg);
	idx->sv_bl = kom_calloc(mb_region_t*, idx->l2b->n_ctg);

	while (fgets(line, sizeof(line), fp)) {
		char *p = line, *ctg, *st_s, *en_s, *endp = 0;
		int tid;
		int64_t st, en, len;
		while (isspace((unsigned char)*p)) ++p;
		if (*p == 0 || *p == '#' || strncmp(p, "track", 5) == 0 || strncmp(p, "browser", 7) == 0)
			continue;
		ctg = strtok(p, " \t\r\n");
		st_s = strtok(0, " \t\r\n");
		en_s = strtok(0, " \t\r\n");
		if (ctg == 0 || st_s == 0 || en_s == 0) {
			++skipped_bad;
			continue;
		}
		errno = 0;
		st = strtoll(st_s, &endp, 10);
		if (errno || endp == st_s) {
			++skipped_bad;
			continue;
		}
		errno = 0;
		en = strtoll(en_s, &endp, 10);
		if (errno || endp == en_s || en <= st) {
			++skipped_bad;
			continue;
		}
		tid = find_ctg(idx->l2b, ctg);
		if (tid < 0) {
			++skipped_ctg;
			continue;
		}
		len = idx->l2b->ctg[tid].len;
		if (st < 0) st = 0;
		if (en > len) en = len;
		if (st >= en) continue;
		append_region(idx, tid, st, en);
		++loaded;
	}
	fclose(fp);
	merge_regions(idx);
	if (kom_verbose >= 2) {
		fprintf(stderr, "[M::%s] loaded %ld HMF SV blacklist intervals from %s", __func__, (long)idx->n_sv_blacklist, fn);
		if (skipped_ctg || skipped_bad)
			fprintf(stderr, " (%ld unmatched contigs, %ld malformed rows skipped)", (long)skipped_ctg, (long)skipped_bad);
		fputc('\n', stderr);
	}
	return loaded;
}

static int sv_blacklist_overlap(const mb_idx_t *idx, const mb_hit_t *r)
{
	int32_t lo, hi, n;
	const mb_region_t *a;
	if (idx == 0 || idx->n_sv_blacklist == 0 || r == 0) return 0;
	if (r->tid < 0 || r->tid >= idx->l2b->n_ctg || r->ts >= r->te) return 0;
	n = idx->sv_bl_n[r->tid];
	if (n == 0) return 0;
	a = idx->sv_bl[r->tid];
	lo = 0, hi = n;
	while (lo < hi) {
		int32_t mid = lo + ((hi - lo) >> 1);
		if (a[mid].en <= r->ts) lo = mid + 1;
		else hi = mid;
	}
	return lo < n && a[lo].st < r->te;
}

void mb_apply_sv_blacklist(const mb_idx_t *idx, const mb_opt_t *opt, int n_regs, mb_hit_t *regs)
{
	int32_t i, cap, n_supp = 0;
	if (idx == 0 || opt == 0 || !(opt->flag & MB_F_HMF_SV_BLACKLIST) || idx->n_sv_blacklist == 0) return;
	cap = opt->sv_blacklist_mapq;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].parent == regs[i].id)
			++n_supp;
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		if (!(r->split || r->split_inv || r->inv || (n_supp > 1 && r->parent == r->id))) continue;
		if (!sv_blacklist_overlap(idx, r)) continue;
		r->sv_blacklist = 1;
		if (cap >= 0 && r->mapq > cap) r->mapq = cap;
	}
}
