#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "kommon.h"
#include "ksort.h"
#include "human_resources.h"

#define key_128x(a) ((a).x)
KRADIX_SORT_INIT(mb128x, mb128_t, key_128x, 8)

#define key_64(a) (a)
KRADIX_SORT_INIT(mb64, uint64_t, key_64, 8)

typedef struct {
	const char *ctg;
	uint64_t st, en;
} mb_builtin_intv_t;

static const mb_builtin_intv_t grch38_problematic_intv[] = {
	{ "chr8",  30393636,  30407636 },
	{ "chr9",  70701717,  70737717 },
	{ "chr11",  1945622,   2041622 },
	{ "chr13", 86253328,  86269328 },
	{ "chr13",111669328, 111703328 },
	{ "chr13",111754328, 111793328 },
	{ "chr16", 19786345,  19801345 },
	{ "chr16", 34339345,  34500345 },
	{ "chr16", 34867345,  34896345 },
	{ "chr16", 34960345,  35072345 },
	{ "chr21",  5009983,   5165983 },
	{ "chr21",  5966983,   6160983 },
	{ "chr21",  6426983,   6579983 },
	{ "chr21",  6789983,   6933983 },
	{ "chr21",  7743983,   7865983 },
	{ "chr21", 13713983,  13718983 },
	{ "chr21", 13752983,  13798983 },
	{ "chr21", 34374983,  34494983 },
	{ "chr21", 43035983,  43186983 },
	{ "chr21", 43376983,  43570983 },
	{ "chr21", 44095983,  44252983 },
	{ "chr22", 18885468,  18939468 },
	{ "chrX",  37084895,  37098895 }
};

static int intv_cmp(const void *a, const void *b)
{
	const l2b_intv_t *x = (const l2b_intv_t*)a, *y = (const l2b_intv_t*)b;
	if (x->st < y->st) return -1;
	if (x->st > y->st) return 1;
	if (x->en < y->en) return -1;
	if (x->en > y->en) return 1;
	return 0;
}

static void mb_idx_clear_problematic(mb_idx_t *idx)
{
	free(idx->problematic);
	idx->problematic = 0;
	idx->n_problematic = idx->m_problematic = 0;
}

static int mb_idx_name2tid(const l2b_t *l2b, const char *name)
{
	uint64_t i;
	for (i = 0; i < l2b->n_ctg; ++i)
		if (strcmp(l2b->ctg[i].name, name) == 0)
			return (int)i;
	if (strncmp(name, "chr", 3) == 0) {
		const char *no_chr = name + 3;
		for (i = 0; i < l2b->n_ctg; ++i)
			if (strcmp(l2b->ctg[i].name, no_chr) == 0)
				return (int)i;
	} else {
		for (i = 0; i < l2b->n_ctg; ++i) {
			const char *ctg = l2b->ctg[i].name;
			if (strncmp(ctg, "chr", 3) == 0 && strcmp(ctg + 3, name) == 0)
				return (int)i;
		}
	}
	return -1;
}

static void mb_idx_add_problematic(mb_idx_t *idx, int32_t tid, uint64_t st, uint64_t en)
{
	const l2b_ctg_t *ctg;
	if (tid < 0 || tid >= (int64_t)idx->l2b->n_ctg) return;
	ctg = &idx->l2b->ctg[tid];
	if (st >= ctg->len) return;
	if (en > ctg->len) en = ctg->len;
	if (st >= en) return;
	kom_grow(l2b_intv_t, idx->problematic, idx->n_problematic, idx->m_problematic);
	idx->problematic[idx->n_problematic].st = ctg->off + st;
	idx->problematic[idx->n_problematic].en = ctg->off + en;
	++idx->n_problematic;
}

static void mb_idx_finish_problematic(mb_idx_t *idx)
{
	uint64_t i, n;
	if (idx->n_problematic == 0) return;
	qsort(idx->problematic, idx->n_problematic, sizeof(*idx->problematic), intv_cmp);
	for (i = 1, n = 0; i < idx->n_problematic; ++i) {
		l2b_intv_t *last = &idx->problematic[n];
		if (idx->problematic[i].st <= last->en) {
			if (last->en < idx->problematic[i].en)
				last->en = idx->problematic[i].en;
		} else {
			++n;
			if (n != i) idx->problematic[n] = idx->problematic[i];
		}
	}
	idx->n_problematic = n + 1;
}

/*****************
 * Index loading *
 *****************/

mb_idx_t *mb_idx_load(const char *prefix, int32_t is_meth)
{
	char *buf;
	mb_idx_t *idx = 0;
	l2b_t *l2b;
	mb_bwt_t *bwt;
	buf = kom_calloc(char, strlen(prefix) + 10);
	strcat(strcpy(buf, prefix), ".l2b");
	l2b = l2b_load(buf);
	if (l2b == 0) goto end_idx_load;
	mb_par_init(l2b);
	if (is_meth) strcat(strcpy(buf, prefix), ".meth.mbw");
	else strcat(strcpy(buf, prefix), ".mbw");
	bwt = mb_bwt_load(buf);
	if (bwt == 0) {
		l2b_destroy(l2b);
		goto end_idx_load;
	}
	mb_bwt_cache(bwt, 10); // TODO: don't hard code this
	mb_human_resources_warn_mismatch(prefix, l2b);
	idx = kom_calloc(mb_idx_t, 1);
	idx->is_meth = !!is_meth, idx->l2b = l2b, idx->bwt = bwt;
end_idx_load:
	free(buf);
	return idx;
}

void mb_idx_destroy(mb_idx_t *idx)
{
	if (idx == 0) return;
	mb_idx_clear_sv_blacklist(idx);
	mb_bwt_destroy(idx->bwt);
	l2b_destroy(idx->l2b);
	free(idx->problematic);
	free(idx);
}

int mb_idx_set_grch38_problematic(mb_idx_t *idx)
{
	uint64_t i;
	if (idx == 0 || idx->l2b == 0) return -1;
	mb_idx_clear_problematic(idx);
	for (i = 0; i < sizeof(grch38_problematic_intv) / sizeof(grch38_problematic_intv[0]); ++i) {
		const mb_builtin_intv_t *r = &grch38_problematic_intv[i];
		int32_t tid = mb_idx_name2tid(idx->l2b, r->ctg);
		if (tid >= 0) mb_idx_add_problematic(idx, tid, r->st, r->en);
	}
	mb_idx_finish_problematic(idx);
	return idx->n_problematic > 0? 0 : -1;
}

int mb_idx_load_problematic_bed(mb_idx_t *idx, const char *fn)
{
	FILE *fp;
	char line[4096];
	if (idx == 0 || idx->l2b == 0 || fn == 0) return -1;
	fp = strcmp(fn, "-") == 0? stdin : fopen(fn, "r");
	if (fp == 0) return -1;
	mb_idx_clear_problematic(idx);
	while (fgets(line, sizeof(line), fp) != 0) {
		char *ctg, *st_s, *en_s, *p = line;
		uint64_t st, en;
		int32_t tid;
		while (*p == ' ' || *p == '\t') ++p;
		if (*p == 0 || *p == '\n' || *p == '#' || strncmp(p, "track", 5) == 0 || strncmp(p, "browser", 7) == 0)
			continue;
		ctg = strtok(p, " \t\r\n");
		st_s = strtok(0, " \t\r\n");
		en_s = strtok(0, " \t\r\n");
		if (ctg == 0 || st_s == 0 || en_s == 0) continue;
		st = strtoull(st_s, 0, 10);
		en = strtoull(en_s, 0, 10);
		if (st >= en) continue;
		tid = mb_idx_name2tid(idx->l2b, ctg);
		if (tid >= 0) mb_idx_add_problematic(idx, tid, st, en);
	}
	if (fp != stdin) fclose(fp);
	mb_idx_finish_problematic(idx);
	return idx->n_problematic > 0? 0 : -1;
}

static int mb_hit_overlap_problematic(const mb_idx_t *idx, const mb_hit_t *r)
{
	uint64_t st, en, lo, hi;
	const l2b_intv_t *a;
	if (idx == 0 || idx->n_problematic == 0 || r == 0 || r->tid < 0 || r->tid >= (int64_t)idx->l2b->n_ctg)
		return 0;
	st = idx->l2b->ctg[r->tid].off + r->ts;
	en = idx->l2b->ctg[r->tid].off + r->te;
	lo = 0, hi = idx->n_problematic;
	while (lo < hi) {
		uint64_t mid = (lo + hi) >> 1;
		if (idx->problematic[mid].en <= st) lo = mid + 1;
		else hi = mid;
	}
	a = lo < idx->n_problematic? &idx->problematic[lo] : 0;
	return a && a->st < en && st < a->en;
}

void mb_apply_problematic_mask(const mb_idx_t *idx, int32_t n_regs, mb_hit_t *regs, int32_t mapq_cap)
{
	int32_t i;
	if (idx == 0 || idx->n_problematic == 0 || regs == 0) return;
	if (mapq_cap < 0) return;
	if (mapq_cap > 60) mapq_cap = 60;
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		if (mb_hit_overlap_problematic(idx, r)) {
			r->problematic = 1;
			if (r->mapq > mapq_cap) r->mapq = mapq_cap;
		}
	}
}

const char *mb_idx_ctg_name(const mb_idx_t *idx, int32_t tid)
{
	return tid >= 0 && tid < idx->l2b->n_ctg? idx->l2b->ctg[tid].name : 0;
}

int64_t mb_idx_ctg_len(const mb_idx_t *idx, int32_t tid)
{
	return tid >= 0 && tid < idx->l2b->n_ctg? idx->l2b->ctg[tid].len : -1;
}

/*****************
 * Thread buffer *
 *****************/

struct mb_tbuf_s {
	void *km;
};

mb_tbuf_t *mb_tbuf_init(int no_kalloc)
{
	mb_tbuf_t *b;
	b = kom_calloc(mb_tbuf_t, 1);
	if (!no_kalloc) b->km = km_init();
	return b;
}

void *mb_tbuf_km(mb_tbuf_t *b)
{
	return b->km;
}

void mb_tbuf_destroy(mb_tbuf_t *b)
{
	if (b->km) km_destroy(b->km);
	free(b);
}

int32_t mb_tbuf_reset(mb_tbuf_t *b, int64_t max_blk_sz)
{
	km_stat_t kmst;
	int64_t max_sz = max_blk_sz < 1U<<28? max_blk_sz : 1U<<28;
	if (b->km == 0) return 0;
	km_stat(b->km, &kmst);
	assert(kmst.n_blocks == kmst.n_cores);
	if (kmst.largest > max_sz || kmst.capacity > max_sz * 2) {
		km_destroy(b->km);
		b->km = km_init();
		return 1;
	}
	return 0;
}

/************************
 * Basic hit operations *
 ************************/

int32_t mb_cal_mblen(int32_t n, const mb_anchor_t *a, int32_t *blen_)
{
	int32_t i;
	int64_t mlen, blen;
	*blen_ = 0;
	if (n <= 0) return 0;
	mlen = blen = a[0].len;
	for (i = 1; i < n; ++i) {
		int span = a[i].len;
		int tl = (int32_t)a[i].tpos - (int32_t)a[i-1].tpos;
		int ql = (int32_t)a[i].qpos - (int32_t)a[i-1].qpos;
		blen += tl > ql? tl : ql;
		mlen += tl > span && ql > span? span : tl < ql? tl : ql;
	}
	*blen_ = blen;
	return mlen;
}

static void mb_cal_fuzzy_len(mb_hit_t *r, const mb_anchor_t *a)
{
	r->mlen = mb_cal_mblen(r->cnt, &a[r->as], &r->blen);
}

static inline void mb_hit_set_coor(mb_hit_t *r, int32_t qlen, const l2b_t *l2b, const mb_anchor_t *a)
{ // NB: r->as and r->cnt MUST BE set correctly for this function to work
	int32_t k = r->as;
	const mb_anchor_t *ak0 = &a[k];
	const mb_anchor_t *ak1 = &a[k + r->cnt - 1];

	r->tid = ak0->sid>>1, r->rev = ak0->sid&1;
	r->ts = ak0->tpos + 1 - ak0->len;
	r->te = ak1->tpos + 1;
	if (!r->rev) { // forward strand
		r->qs = ak0->qpos + 1 - ak0->len;
		r->qe = ak1->qpos + 1;
	} else { // reverse strand
		r->qs = qlen - (ak1->qpos + 1);
		r->qe = qlen - (ak0->qpos + 1 - ak0->len);
	}
	mb_cal_fuzzy_len(r, a);
}

int32_t mb_cal_high_cov(void *km, int32_t n, const mb_sai_t *sai, int32_t max_occ)
{
	int32_t i, n_hi = 0, hi_st, hi_en, hi_cov;
	uint64_t *b;
	for (i = 0; i < n; ++i)
		if (sai[i].size > max_occ)
			++n_hi;
	if (n_hi == 0) return 0;
	b = Kmalloc(km, uint64_t, n_hi);
	for (i = 0, n_hi = 0; i < n; ++i)
		if (sai[i].size > max_occ)
			b[n_hi++] = sai[i].info;
	radix_sort_mb64(b, b + n_hi);
	hi_st = b[0]>>32, hi_en = (int32_t)b[0], hi_cov = 0;
	for (i = 1; i < n_hi; ++i) {
		int32_t st = b[i]>>32, en = (int32_t)b[i];
		if (st > hi_en) {
			hi_cov += hi_en - hi_st;
			hi_st = st, hi_en = en;
		} else hi_en = hi_en > en? hi_en : en;
	}
	hi_cov += hi_en - hi_st;
	kfree(km, b);
	return hi_cov;
}

void mb_sync_high_cov(int32_t n, mb_hit_t *h)
{
	int32_t i, max_frac = 0;
	for (i = 0; i < n; ++i)
		max_frac = max_frac > h[i].frac_high? max_frac : h[i].frac_high;
	for (i = 0; i < n; ++i)
		h[i].frac_high = max_frac;
}

mb_hit_t *mb_gen_hit(void *km, uint32_t hash, int qlen, const l2b_t *l2b, int n_u, uint64_t *u, mb_anchor_t *a)
{ // convert chains to hits
	mb128_t *z, tmp;
	mb_hit_t *r;
	int i, k;

	if (n_u <= 0) return 0;

	// sort by score
	z = Kmalloc(km, mb128_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		uint32_t h;
		h = (uint32_t)mb_hash64((mb_hash64(a[k].tpos) + mb_hash64(a[k].qpos)) ^ hash);
		z[i].x = u[i] ^ h; // u[i] -- higher 32 bits: chain score; lower 32 bits: number of anchors
		z[i].y = (uint64_t)k << 32 | (int32_t)u[i];
		k += (int32_t)u[i];
	}
	radix_sort_mb128x(z, z + n_u);
	for (i = 0; i < n_u>>1; ++i) // reverse, s.t. larger score first
		tmp = z[i], z[i] = z[n_u-1-i], z[n_u-1-i] = tmp;

	// populate r[]
	r = (mb_hit_t*)calloc(n_u, sizeof(mb_hit_t));
	for (i = 0; i < n_u; ++i) {
		mb_hit_t *ri = &r[i];
		ri->id = i;
		ri->parent = MB_PARENT_UNSET;
		ri->score = ri->score0 = z[i].x >> 32;
		ri->hash = (uint32_t)z[i].x;
		ri->cnt = (int32_t)z[i].y;
		ri->as = z[i].y >> 32;
		mb_hit_set_coor(ri, qlen, l2b, a);
	}
	kfree(km, z);
	return r;
}

void mb_split_hit(mb_hit_t *r, mb_hit_t *r2, int n, int qlen, mb_anchor_t *a, const l2b_t *l2b)
{
	if (n <= 0 || n >= r->cnt) return;
	*r2 = *r;
	r2->id = -1;
	r2->sam_pri = 0;
	r2->p = 0;
	r2->split_inv = 0;
	r2->cnt = r->cnt - n;
	r2->score = (int32_t)(r->score * ((float)r2->cnt / r->cnt) + .499);
	r2->as = r->as + n;
	if (r->parent == r->id) r2->parent = MB_PARENT_TMP_PRI;
	mb_hit_set_coor(r2, qlen, l2b, a);
	r->cnt -= r2->cnt;
	r->score -= r2->score;
	mb_hit_set_coor(r, qlen, l2b, a);
	r->split |= 1, r2->split |= 2;
}

void mb_sync_hits(void *km, int n_regs, mb_hit_t *regs)
{
	int *tmp, i, max_id = -1, n_tmp;
	if (n_regs <= 0) return;
	for (i = 0; i < n_regs; ++i)
		max_id = max_id > regs[i].id? max_id : regs[i].id;
	n_tmp = max_id + 1;
	tmp = (int*)kmalloc(km, n_tmp * sizeof(int));
	for (i = 0; i < n_tmp; ++i) tmp[i] = -1;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].id >= 0) tmp[regs[i].id] = i;
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		r->id = i;
		if (r->parent == MB_PARENT_TMP_PRI)
			r->parent = i;
		else if (r->parent >= 0 && tmp[r->parent] >= 0)
			r->parent = tmp[r->parent];
		else r->parent = MB_PARENT_UNSET;
	}
	kfree(km, tmp);
	mb_set_sam_pri(n_regs, regs, 0); // this flag will be overwritten later anyway
}

/**********************************
 * Set primary and secondary hits *
 **********************************/

static int update_sub(const l2b_t *l2b, mb_hit_t *ri, mb_hit_t *rp, float mask_level, int mask_len, int sub_diff, int uncov_len)
{
	int si = ri->qs, ei = ri->qe, sj = rp->qs, ej = rp->qe, min, max, ol;
	if (ej <= si || sj >= ei) return 0;
	min = ej - sj < ei - si? ej - sj : ei - si;
	max = ej - sj > ei - si? ej - sj : ei - si;
	ol = ej <= si || sj >= ei? 0 : (ej < ei? ej : ei) - (sj > si? sj : si);
	if ((double)ol / min - (double)uncov_len / max > mask_level && uncov_len <= mask_len) {
		int cnt_sub = 0, sci = ri->score;
		ri->parent = rp->parent;
		if (mb_par_equiv(l2b, ri, rp)) return 1;
		rp->subsc = rp->subsc > sci? rp->subsc : sci;
		if (rp->p && ri->p && (rp->tid != ri->tid || rp->ts != ri->ts || rp->te != ri->te || ol != min)) { // the last condition excludes identical hits after DP
			sci = ri->p->dp_max;
			rp->p->dp_max2 = rp->p->dp_max2 > sci? rp->p->dp_max2 : sci;
			if (rp->p->dp_max - ri->p->dp_max <= sub_diff) cnt_sub = 1;
		}
		if (cnt_sub) ++rp->n_sub;
		return 1;
	} else return 0;
}

void mb_set_parent(void *km, const l2b_t *l2b, float mask_level, int mask_len, int n, mb_hit_t *r, int sub_diff, int hard_mask_level)
{ // TODO: re-examine the logic for variable-length seeds
	int i, j, k, *w;
	uint64_t *cov;
	if (n <= 0) return;
	for (i = 0; i < n; ++i) r[i].id = i;
	cov = Kmalloc(km, uint64_t, n);
	w = Kmalloc(km, int, n);
	w[0] = 0, r[0].parent = 0;
	for (i = 1, k = 1; i < n; ++i) {
		mb_hit_t *ri = &r[i];
		int si = ri->qs, ei = ri->qe, n_cov = 0, uncov_len = 0, max_ol, max_j, n_par = 0;
		if (hard_mask_level) goto skip_uncov;
		for (j = 0; j < k; ++j) {
			const mb_hit_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe;
			if (ej <= si || sj >= ei) continue;
			if (sj < si) sj = si;
			if (ej > ei) ej = ei;
			cov[n_cov++] = (uint64_t)sj<<32 | ej;
		}
		if (n_cov == 0) {
			goto add_primary;
		} else {
			int j, x = si;
			radix_sort_mb64(cov, cov + n_cov);
			for (j = 0; j < n_cov; ++j) {
				if ((int)(cov[j]>>32) > x) uncov_len += (cov[j]>>32) - x;
				x = (int32_t)cov[j] > x? (int32_t)cov[j] : x;
			}
			if (ei > x) uncov_len += ei - x;
		}
skip_uncov:
		for (j = 0, max_ol = 0, max_j = -1; j < k; ++j) { // find the parent with the maximum overlap
			const mb_hit_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe;
			int ol = ej <= si || sj >= ei? 0 : (ej < ei? ej : ei) - (sj > si? sj : si);
			if (max_ol < ol) max_ol = ol, max_j = j;
		}
		if (max_j >= 0) {
			n_par += update_sub(l2b, ri, &r[w[max_j]], mask_level, mask_len, sub_diff, uncov_len);
			for (j = 0; j < k && n_par == 0; ++j) // if no parent found on the longest overlap, try more
				n_par += update_sub(l2b, ri, &r[w[j]], mask_level, mask_len, sub_diff, uncov_len);
		}
add_primary:
		if (n_par == 0) w[k++] = i, ri->parent = i, ri->n_sub = 0;
	}
	kfree(km, cov);
	kfree(km, w);
}

void mb_set_sam_pri(int32_t n, mb_hit_t *r, int32_t is_primary5)
{
	int32_t i, n_pri = 0, min_i = -1, min_qs = -1, first_i = -1;
	if (n <= 0) return;
	for (i = 0; i < n; ++i) {
		r[i].sam_pri = 0;
		if (r[i].id != r[i].parent) continue;
		if (++n_pri == 1) first_i = i;
		if (min_qs < 0 || r[i].qs < min_qs)
			min_i = i, min_qs = r[i].qs;
	}
	assert(n_pri > 0);
	if (is_primary5) r[min_i].sam_pri = 1;
	else r[first_i].sam_pri = 1;
}

void mb_select_sub(void *km, float pri_ratio, int min_diff, int best_n, int *n_, mb_hit_t *r)
{
	if (pri_ratio > 0.0f && *n_ > 0) {
		int i, k, n = *n_, n_2nd = 0;
		uint8_t *keep = Kcalloc(km, uint8_t, n);
		for (i = 0; i < n; ++i) {
			int p = r[i].parent;
			if (p == i || r[i].inv) {
				keep[i] = 1;
			} else if ((r[i].score >= r[p].score * pri_ratio || r[i].score + min_diff >= r[p].score) && n_2nd < best_n) {
				if (!(r[i].qs == r[p].qs && r[i].qe == r[p].qe && r[i].tid == r[p].tid && r[i].ts == r[p].ts && r[i].te == r[p].te))
					keep[i] = 1, ++n_2nd;
			}
		}
		for (i = k = 0; i < n; ++i) {
			if (keep[i]) r[k++] = r[i];
			else if (r[i].p) free(r[i].p); // r->p is libc-allocated; free here before the pointer is lost
		}
		kfree(km, keep);
		if (k != n) mb_sync_hits(km, k, r);
		*n_ = k;
	}
}

void mb_hit_sort(void *km, int *n_regs, mb_hit_t *r)
{
	int32_t i, n_aux, n = *n_regs;
	mb128_t *aux;
	mb_hit_t *t;

	if (n <= 1) return;
	aux = (mb128_t*)kmalloc(km, (size_t)n * 16);
	t = (mb_hit_t*)kmalloc(km, (size_t)n * sizeof(mb_hit_t));
	for (i = n_aux = 0; i < n; ++i) {
		if (r[i].inv || r[i].cnt >= 0) {
			int score = r[i].p? r[i].p->dp_max : r[i].score;
			aux[n_aux].x = (uint64_t)score << 32 | r[i].hash;
			aux[n_aux++].y = i;
		} else if (r[i].p) {
			free(r[i].p);
			r[i].p = 0;
		}
	}
	radix_sort_mb128x(aux, aux + n_aux);
	for (i = n_aux - 1; i >= 0; --i)
		t[n_aux - 1 - i] = r[aux[i].y];
	memcpy(r, t, sizeof(mb_hit_t) * n_aux);
	*n_regs = n_aux;
	kfree(km, aux);
	kfree(km, t);
}

void mb_filter_hits(const mb_opt_t *opt, int qlen, int *n_regs, mb_hit_t *regs)
{
	int i, k;
	for (i = k = 0; i < *n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		int flt = r->flt;
		if (r->p) {
			if (r->mlen < opt->min_chain_score) flt = 1;
			else if (r->p->dp_max < opt->min_dp_max * opt->a) flt = 1;
			if (flt) free(r->p);
		}
		if (!flt) {
			if (k < i) regs[k++] = regs[i];
			else ++k;
		}
	}
	*n_regs = k;
}

int mb_squeeze_a(void *km, int n_regs, mb_hit_t *regs, mb_anchor_t *a)
{
	int i, as = 0;
	uint64_t *aux;
	aux = (uint64_t*)kmalloc(km, (size_t)n_regs * 8);
	for (i = 0; i < n_regs; ++i)
		aux[i] = (uint64_t)regs[i].as << 32 | i;
	radix_sort_mb64(aux, aux + n_regs);
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[(int32_t)aux[i]];
		if (r->as != as) {
			memmove(&a[as], &a[r->as], (size_t)r->cnt * sizeof(mb_anchor_t));
			r->as = as;
		}
		as += r->cnt;
	}
	kfree(km, aux);
	return as;
}

/*******************
 * Mapping quality *
 *******************/

static void mb_set_inv_mapq(void *km, int n_regs, mb_hit_t *regs)
{
	int i, n_aux;
	mb128_t *aux;
	if (n_regs < 3) return;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].inv) break;
	if (i == n_regs) return; // no inversion hits

	aux = Kmalloc(km, mb128_t, n_regs);
	for (i = n_aux = 0; i < n_regs; ++i)
		if (regs[i].parent == i || regs[i].parent < 0)
			aux[n_aux].y = i, aux[n_aux++].x = (uint64_t)regs[i].tid << 32 | regs[i].ts;
	radix_sort_mb128x(aux, aux + n_aux);

	for (i = 1; i < n_aux - 1; ++i) {
		mb_hit_t *inv = &regs[aux[i].y];
		if (inv->inv) {
			mb_hit_t *l = &regs[aux[i-1].y];
			mb_hit_t *r = &regs[aux[i+1].y];
			inv->mapq = l->mapq < r->mapq? l->mapq : r->mapq;
		}
	}
	kfree(km, aux);
}

void mb_cap_mapq_by_mask(const l2b_t *l2b, int n_regs, mb_hit_t *regs, float mask_level)
{
	int i;
	if (l2b == 0 || l2b->n_mask == 0 || mask_level <= 0.0f) return;
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		int64_t len, masked;
		if (r->mapq == 0 || r->tid < 0 || r->te <= r->ts) continue;
		len = r->te - r->ts;
		masked = l2b_mask_overlap(l2b, r->tid, r->ts, r->te);
		if (masked > mask_level * len && r->mapq > 10)
			r->mapq = 10;
	}
}

void mb_set_mapq(void *km, const l2b_t *l2b, int32_t qlen, int n_regs, mb_hit_t *regs, int min_chain_sc, int match_sc, int is_sr, int max_sr_len, float mask_level)
{
	const int32_t mapQ_coef_len = 50;
	const double mapQ_coef_fac = 3.0; // should be log(mapQ_coef_len)), but bwa-mem uses 3.0 due to a bug. Let's match bwa-mem
	const double q_coef = 40.0f;
	int i;
	if (n_regs == 0) return;
	for (i = 0; i < n_regs; ++i) {
		mb_hit_t *r = &regs[i];
		if (r->inv) {
			r->mapq = 0;
		} else if (r->parent == r->id) {
			int mapq, mapq_sr, mapq_lr, subsc;
			double pen_chn = r->score > qlen * 0.1? 1.0 : 10.0 * r->score / qlen; // penalize chains with few matching bases
			subsc = r->subsc > min_chain_sc? r->subsc : min_chain_sc;
			if (r->p && r->p->dp_max2 > 0 && r->p->dp_max > 0) {
				double x, identity = (double)r->mlen / r->blen;
				// BWA-MEM formula for short reads
				x = r->blen < mapQ_coef_len? 1. : mapQ_coef_fac / log(r->blen);
				x *= identity * identity;
				mapq_sr = (int)(6.02 * x * x * (r->p->dp_max - r->p->dp_max2) / match_sc + .499f);
				// minimap2 formula for long reads
				x = (double)r->p->dp_max2 / r->p->dp_max;
				if (subsc > r->score0) x *= (double)subsc / r->score0;
				mapq_lr = (int)(pen_chn * identity * q_coef * (1.0 - x * x) * log((double)r->p->dp_max / match_sc));
				// final mapq
				if (is_sr) mapq = mapq_sr;
				else if (max_sr_len < 0) mapq = mapq_lr;
				else mapq = qlen < max_sr_len? mapq_sr : (int32_t)(mapq_lr - (mapq_lr - mapq_sr) * pow(2.0, 1.0 - (double)qlen / max_sr_len) + .499);
			} else { // minimap2 formula
				double x = (double)subsc / r->score0;
				if (r->p) {
					double identity = (double)r->mlen / r->blen;
					mapq = (int)(pen_chn * identity * q_coef * (1.0f - x) * log((double)r->p->dp_max / match_sc));
				} else {
					mapq = (int)(pen_chn * q_coef * (1.0f - x) * log(r->score));
				}
			}
			mapq -= (int)(4.343f * log(r->n_sub + 1) + .499f);
			mapq = mapq > 0? mapq : 0;
			if (r->seed_ratio < 50) mapq *= (double)r->seed_ratio * r->seed_ratio / 2500.0;
			r->mapq = mapq < 60? mapq : 60;
			if (r->p && r->p->dp_max > r->p->dp_max2 && r->mapq == 0) r->mapq = 1;
		} else r->mapq = 0;
	}
	mb_cap_mapq_by_mask(l2b, n_regs, regs, mask_level);
	mb_set_inv_mapq(km, n_regs, regs);
}

/************************
 * Core mapping routine *
 ************************/

static void mb_dbg_seed(int64_t n, const mb_sai_t *u, const char *qname)
{
	int64_t i;
	for (i = 0; i < n; ++i) {
		const mb_sai_t *p = &u[i];
		fprintf(stderr, "SD\t%s\t%d\t%d\t%ld\n", qname? qname : "*", (int32_t)(p->info>>32), (int32_t)p->info, (long)p->size);
	}
}

static void mb_dbg_anchor(const mb_idx_t *idx, int qlen, int64_t n, const mb_anchor_t *a, const char *qname)
{
	int64_t i;
	for (i = 0; i < n; ++i) {
		const mb_anchor_t *ai = &a[i];
		int rid = ai->sid >> 1;
		int rev = ai->sid & 1;
		int32_t qs = rev? qlen - 1 - ai->qpos : ai->qpos + 1 - ai->len;
		int64_t ts = ai->tpos + 1 - ai->len;
		fprintf(stderr, "AC\t%s\t%d\t%c\t%s\t%ld\t%d\n", qname? qname : "*", qs, "+-"[rev], idx->l2b->ctg[rid].name, (long)ts, ai->len);
	}
}

static mb_hit_t *mb_map_sai_core(const mb_opt_t *opt, const mb_idx_t *idx, int64_t qlen, const char *seq0, const uint8_t *seq4, l2b_meth_t mt, mb_sai_v *u, int32_t *n_hit_, mb_tbuf_t *b, const char *qname)
{
	const int32_t min_rechain_len = 1000;
	const double min_rechain_ratio = 0.1;
	uint32_t hash;
	int32_t i, n_hit, hi_cov, is_sr;
	int32_t sub_diff = opt->a + opt->b > opt->q + opt->e? opt->a + opt->b : opt->q + opt->e;
	uint64_t *w;
	double chn_pen_gap, seed_ratio;
	uint8_t *seq;
	mb_anchor_v v = {0,0,0};
	mb_anchor_t *a;
	mb_hit_t *hit;

	if (kom_dbg_flag & MB_DBG_QNAME) fprintf(stderr, "QN\t%s\n", qname);

	*n_hit_ = 0;
	if (u->n == 0) {
		kfree(b->km, u->a);
		return 0;
	}
	hash  = qname? mb_hash_str(qname) : 0;
	hash ^= mb_hash64(qlen) + mb_hash64(opt->seed);
	hash  = mb_hash64(hash);
	if (seq4) {
		seq = (uint8_t*)seq4;
	} else {
		seq = kmalloc(b->km, qlen);
		for (i = 0; i < qlen; ++i) seq[i] = kom_nt4_table[(uint8_t)seq0[i]];
	}
	hi_cov = mb_cal_high_cov(b->km, u->n, u->a, opt->max_occ);
	is_sr = mb_is_sr_mode(opt, qlen);

	// collect anchors
	chn_pen_gap = opt->chain_gap_scale * .01 * opt->min_len;
	if (kom_dbg_flag & MB_DBG_SEED) mb_dbg_seed(u->n, u->a, qname);
	seed_ratio = mb_anchor(b->km, idx, u, opt->min_len, qlen, seq, mt, opt->max_occ, opt->max_sub_occ, &v);
	kfree(b->km, u->a); // no longer needed
	u->n = 0, u->a = 0;

	// initial chaining
	if (kom_dbg_flag & MB_DBG_ANCHOR) mb_dbg_anchor(idx, qlen, v.n, v.a, qname);
	a = mb_lchain_dp(b->km, idx->l2b, opt->max_gap, opt->max_gap, opt->bw, opt->max_chain_skip, opt->max_chain_iter,
					 opt->min_chain_score, chn_pen_gap, v.n, v.a, &n_hit, &w);
	v.a = 0; v.n = v.m = 0; // ownership transferred to a

	// re-chaining
	if (opt->bw_long > opt->bw * 2 && !is_sr && n_hit > 0) {
		int64_t n_a, as, st, en;
		int32_t best;
		// chains in w[] are sorted by tpos of first anchor, not by score; find the best
		for (i = 1, best = 0; i < n_hit; ++i)
			if ((w[i] >> 32) > (w[best] >> 32)) best = i;
		for (i = 0, as = 0; i < best; ++i) as += (int32_t)w[i];
		st = a[as].qpos + 1 - a[as].len;
		en = a[as + (int32_t)w[best] - 1].qpos + 1;
		if (qlen - (en - st) > min_rechain_len && en - st > qlen * min_rechain_ratio) {
			for (i = 0, n_a = 0; i < n_hit; ++i) n_a += (int32_t)w[i];
			kfree(b->km, w);
			mb_anchor_sort(idx->l2b, n_a, a);
			a = mb_lchain_dp(b->km, idx->l2b, opt->max_gap, opt->max_gap, opt->bw_long, opt->max_chain_skip, opt->max_chain_iter,
							 opt->min_chain_score, chn_pen_gap, n_a, a, &n_hit, &w);
		}
	}

	// chain ordering
	hit = mb_gen_hit(b->km, hash, qlen, idx->l2b, n_hit, w, a);
	kfree(b->km, w);
	mb_mark_par_hits(idx->l2b, n_hit, hit);
	mb_set_parent(b->km, idx->l2b, opt->mask_level, opt->mask_len, n_hit, hit, sub_diff, 0);
	mb_par_resolve(idx->l2b, n_hit, hit, sub_diff);
	mb_select_sub(b->km, opt->pri_ratio, opt->min_len * 2, opt->best_n, &n_hit, hit);

	// base alignment
	if (!(opt->flag & MB_F_NO_ALN)) {
		hit = mb_align_skeleton(b->km, opt, idx, qlen, seq, mt, &n_hit, hit, a);
		mb_mark_par_hits(idx->l2b, n_hit, hit);
		mb_set_parent(b->km, idx->l2b, opt->mask_level, opt->mask_len, n_hit, hit, sub_diff, 0);
		mb_par_resolve(idx->l2b, n_hit, hit, sub_diff);
		mb_select_sub(b->km, opt->pri_ratio, opt->min_len * 2, opt->best_n, &n_hit, hit);
		mb_set_sam_pri(n_hit, hit, !!(opt->flag & MB_F_PRIMARY5));
	}
	for (i = 0; i < n_hit; ++i) {
		hit[i].frac_high = (int32_t)(255. * hi_cov / qlen);
		hit[i].seed_ratio = (int32_t)(255. * seed_ratio + .499);
		if (hit[i].seed_ratio == 0) hit[i].seed_ratio = 1;
	}
	mb_set_mapq(b->km, idx->l2b, qlen, n_hit, hit, opt->min_chain_score, opt->a, is_sr, opt->max_sr_len, opt->mask_level);
	mb_apply_sv_blacklist(idx, opt, n_hit, hit);
	if (opt->flag & MB_F_PROBLEMATIC_MASK)
		mb_apply_problematic_mask(idx, n_hit, hit, opt->problematic_mapq_cap);
	mb_apply_unmap_regions(opt->unmap_regions, n_hit, hit);
	mb_mapq_track_apply(opt->mapq_track, n_hit, hit);

	// clean up
	kfree(b->km, a);
	if (seq4 == 0) kfree(b->km, seq);
	*n_hit_ = n_hit;
	return hit;
}

mb_hit_t *mb_map_sai(const mb_opt_t *opt, const mb_idx_t *idx, int64_t qlen, const char *seq0, l2b_meth_t mt, mb_sai_v *u, int32_t *n_hit_, mb_tbuf_t *b, const char *qname)
{
	return mb_map_sai_core(opt, idx, qlen, seq0, 0, mt, u, n_hit_, b, qname);
}

mb_hit_t *mb_map_sai4(const mb_opt_t *opt, const mb_idx_t *idx, int64_t qlen, const uint8_t *seq, l2b_meth_t mt, mb_sai_v *u, int32_t *n_hit_, mb_tbuf_t *b, const char *qname)
{
	return mb_map_sai_core(opt, idx, qlen, 0, seq, mt, u, n_hit_, b, qname);
}

/*************************
 * Public alignment APIs *
 *************************/

mb_hit_t *mb_map(const mb_opt_t *opt, const mb_idx_t *idx, int32_t qlen, const char *seq0, int32_t mt0, int32_t *n_hit_, mb_tbuf_t *b0, const char *qname)
{
	mb_opt_t opt_adap;
	mb_hit_t *ret;
	mb_sai_v u = {0,0,0};
	mb_tbuf_t *b;
	uint8_t *seq;
	int32_t i;
	l2b_meth_t mt = mt0 == 0? L2B_METH_NONE : mt0 == 1? L2B_METH_C2T : L2B_METH_G2A;
	b = b0? b0 : mb_tbuf_init(1);
	mb_opt_adap(opt, qlen, &opt_adap);
	seq = Kmalloc(b->km, uint8_t, qlen);
	for (i = 0; i < qlen; ++i)
		seq[i] = kom_nt4_table[(uint8_t)seq0[i]];
	mb_seed_intv(b->km, idx->bwt, qlen, seq, opt->min_len, opt->max_sub_occ, &u);
	ret = mt == L2B_METH_NONE? mb_map_sai4(&opt_adap, idx, qlen, seq, mt, &u, n_hit_, b, qname)
							  : mb_map_sai(&opt_adap, idx, qlen, seq0, mt, &u, n_hit_, b, qname);
	kfree(b->km, seq);
	if (b0 == 0) mb_tbuf_destroy(b);
	return ret;
}

mb_hit_t **mb_map_batch(const mb_opt_t *opt, const mb_idx_t *idx, int32_t n_seq, const int32_t *qlen, const char **seq, int32_t *n_hit, mb_tbuf_t *b0, const char **qname)
{
	mb_tbuf_t *b;
	mb_hit_t **hit;
	mb_sai_v *sai;
	uint8_t **seq4;
	void *km;
	int32_t i, j, k, sb_st, sb_len, sb_max, is_pe = !!(opt->flag & MB_F_PE);

	if (n_seq <= 0) return 0;
	b = b0? b0 : mb_tbuf_init(0);
	km = mb_tbuf_km(b);
	hit = (mb_hit_t**)calloc(n_seq, sizeof(mb_hit_t*));

	// pre-allocate for sub-batch processing
	sb_max = opt->sb_seq < n_seq? opt->sb_seq : n_seq;
	seq4 = Kmalloc(km, uint8_t*, sb_max);
	sai = Kmalloc(km, mb_sai_v, sb_max);

	// process in sub-batches
	for (i = 0, sb_st = 0, sb_len = 0; i <= n_seq; ++i) {
		if (i == n_seq || sb_len >= opt->sb_len || i - sb_st >= opt->sb_seq) {
			int32_t sb_n = i - sb_st;
			if (sb_n <= 0) { sb_st = i; sb_len = 0; continue; }

			// convert sub-batch to 4-bit encoding
			for (k = 0; k < sb_n; ++k) {
				int32_t idx_k = sb_st + k;
				seq4[k] = Kmalloc(km, uint8_t, qlen[idx_k]);
				for (j = 0; j < qlen[idx_k]; ++j)
					seq4[k][j] = kom_nt4_table[(uint8_t)seq[idx_k][j]];
			}

			// batch SMEM for sub-batch
			memset(sai, 0, (size_t)sb_n * sizeof(*sai));
			mb_seed_intv_batch(km, idx->bwt, sb_n, &qlen[sb_st], seq4, opt->min_len, opt->max_sub_occ, sai);

			// map each sequence in sub-batch
			for (k = 0; k < sb_n; ++k) {
				int32_t idx_k = sb_st + k;
				mb_opt_t opt_adap;
				l2b_meth_t mt = L2B_METH_NONE;
				mb_opt_adap(opt, qlen[idx_k], &opt_adap);
				if (opt->flag & MB_F_METH) {
					if (is_pe) mt = (idx_k&1) == 0? L2B_METH_C2T : L2B_METH_G2A;
					else mt = L2B_METH_C2T;
				}
				if (mt == L2B_METH_NONE)
					hit[idx_k] = mb_map_sai4(&opt_adap, idx, qlen[idx_k], seq4[k], mt, &sai[k], &n_hit[idx_k], b, qname? qname[idx_k] : 0);
				else
					hit[idx_k] = mb_map_sai(&opt_adap, idx, qlen[idx_k], seq[idx_k], mt, &sai[k], &n_hit[idx_k], b, qname? qname[idx_k] : 0);
			}
			for (k = 0; k < sb_n; ++k) kfree(km, seq4[k]);

			sb_st = i;
			sb_len = 0;
		}
		if (i < n_seq) sb_len += qlen[i];
	}

	kfree(km, sai);
	kfree(km, seq4);

	// paired-end processing
	if (is_pe && n_seq >= 2) {
		mb_pestat_t pes[4];
		for (i = 0; i < 4; ++i) pes[i].failed = 1;
		pes[1].failed = 0;
		pes[1].avg = opt->pe_avg, pes[1].std = opt->pe_std;
		pes[1].lo = opt->pe_lo, pes[1].hi = opt->pe_hi;
		for (i = 0; i + 1 < n_seq; i += 2) {
			int32_t len2[2] = { qlen[i], qlen[i+1] };
			char *seq2[2] = { (char*)seq[i], (char*)seq[i+1] };
			mb_pair(km, opt, idx->l2b, &n_hit[i], &hit[i], pes, len2, seq2);
			mb_apply_sv_blacklist(idx, opt, n_hit[i], hit[i]);
			mb_apply_sv_blacklist(idx, opt, n_hit[i+1], hit[i+1]);
			if (opt->flag & MB_F_PROBLEMATIC_MASK) {
				mb_apply_problematic_mask(idx, n_hit[i], hit[i], opt->problematic_mapq_cap);
				mb_apply_problematic_mask(idx, n_hit[i+1], hit[i+1], opt->problematic_mapq_cap);
			}
		}
	}

	if (b0 == 0) mb_tbuf_destroy(b);
	return hit;
}
