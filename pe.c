#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "ksw2.h"

static inline int mb_insert_dir(const mb_hit_t *h0, const mb_hit_t *h1, int64_t *dist)
{
	int64_t p0, p1;
	p0 = h0->rev? h0->te : h0->ts;
	p1 = h1->rev? h1->te : h1->ts;
	*dist = p0 > p1? p0 - p1 : p1 - p0;
	return ((int32_t)h0->rev << 1 | (int32_t)h1->rev) ^ (p0 < p1? 0 : 3);
}

static const mb_hit_t *mb_select_unique_se(int32_t n_hit, const mb_hit_t *hit)
{
	int32_t j, n_pri = 0, mapq = 0, k = -1;
	for (j = 0; j < n_hit; ++j)
		if (hit[j].id == hit[j].parent) 
			++n_pri, mapq = hit[j].mapq, k = j;
	return n_pri == 1 && mapq >= 10? &hit[k] : 0;
}

void mb_pestat(void *km, const mb_opt_t *opt, int32_t n_frag, const int32_t *seg_off, const int32_t *seg_cnt, const int32_t *n_hit, mb_hit_t *const *hit, mb_pestat_t pes[4])
{
	const int MIN_DIR_CNT = 10;
	const double MIN_DIR_RATIO = 0.05, OUTLIER_BOUND = 2.0, MAPPING_BOUND = 3.0, MAX_STDDEV = 4.0;
	int32_t i, d, max;
	struct { int32_t n, m; uint64_t *a; } is[4], *q;
	memset(is, 0, sizeof(is[0]) * 4);
	memset(pes, 0, sizeof(pes[0]) * 4);
	for (i = 0; i < n_frag; ++i) {
		const mb_hit_t *r[2];
		int32_t off, dir;
		int64_t dist;
		if (seg_cnt[i] != 2) continue;
		off = seg_off[i];
		r[0] = mb_select_unique_se(n_hit[off + 0], hit[off + 0]);
		r[1] = mb_select_unique_se(n_hit[off + 1], hit[off + 1]);
		if (r[0] == 0 || r[1] == 0) continue;
		if (r[0]->tid != r[1]->tid) continue; // not on the same contig
		dir = mb_insert_dir(r[0], r[1], &dist);
		if (dist < opt->max_pe_ins) {
			if (is[dir].n == is[dir].m)
				Kgrow(km, uint64_t, is[dir].a, is[dir].n, is[dir].m);
			is[dir].a[is[dir].n++] = dist;
		}
	}
	if (kom_verbose >= 3)
		fprintf(stderr, "[M::%s] # candidate unique pairs for (FF, FR, RF, RR): (%d, %d, %d, %d)\n", __func__, is[0].n, is[1].n, is[2].n, is[3].n);
	for (d = 0, max = 0; d < 4; ++d)
		max = max > is[d].n? max : is[d].n;
	for (d = 0; d < 4; ++d) {
		mb_pestat_t *r = &pes[d];
		q = &is[d];
		int p25, p50, p75, x;
		if (q->n < MIN_DIR_CNT || q->n < max * MIN_DIR_RATIO) {
			r->failed = 1;
			kfree(km, q->a);
			continue;
		}
		radix_sort_mb64(q->a, q->a + q->n);
		p25 = q->a[(int)(.25 * q->n + .499)];
		p50 = q->a[(int)(.50 * q->n + .499)];
		p75 = q->a[(int)(.75 * q->n + .499)];
		r->lo  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
		if (r->lo < 1) r->lo = 1;
		r->hi = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
		for (i = x = 0, r->avg = 0; i < q->n; ++i)
			if (q->a[i] >= r->lo && q->a[i] <= r->hi)
				r->avg += q->a[i], ++x;
		r->avg /= x;
		for (i = 0, r->std = 0; i < q->n; ++i)
			if (q->a[i] >= r->lo && q->a[i] <= r->hi)
				r->std += (q->a[i] - r->avg) * (q->a[i] - r->avg);
		r->std = sqrt(r->std / x);
		if (kom_verbose >= 3)
			fprintf(stderr, "[M::%s::%c%c] (25, 50, 75) percentile: (%d, %d, %d); mean and std.dev: (%.2f, %.2f)\n",
				__func__, "FR"[d>>1&1], "FR"[d&1], p25, p50, p75, r->avg, r->std);
		r->lo  = (int)(p25 - MAPPING_BOUND * (p75 - p25) + .499);
		r->hi = (int)(p75 + MAPPING_BOUND * (p75 - p25) + .499);
		if (r->lo  > r->avg - MAX_STDDEV * r->std) r->lo  = (int)(r->avg - MAX_STDDEV * r->std + .499);
		if (r->hi < r->avg + MAX_STDDEV * r->std) r->hi = (int)(r->avg + MAX_STDDEV * r->std + .499);
		if (r->lo < 1) r->lo = 1;
		if (kom_verbose >= 3)
			fprintf(stderr, "[M::%s::%c%c] low and high boundaries for proper pairs: (%d, %d)\n", __func__, "FR"[d>>1&1], "FR"[d&1], r->lo, r->hi);
		kfree(km, q->a);
	}
}

#define MB_SQRT1_2 0.707106781186547524401

typedef struct {
	int32_t score, sub_sc, n_sub, n_pp;
	int32_t i[2];
} mb_pairaux_t;

static void mb_pair_hits(void *km, const mb_opt_t *opt, const l2b_t *l2b, int32_t n_hit[2], mb_hit_t *hit[2], const mb_pestat_t pes[4], mb_pairaux_t *ret)
{
	int32_t r, i, k, n_pa, y[4], n_pp = 0, m_pp = 0;
	mb128_t *pa, *pp = 0; // pp: proper pairs

	ret->i[0] = ret->i[1] = ret->score = ret->sub_sc = -1, ret->n_sub = 0;
	if (n_hit[0] == 0 || n_hit[1] == 0) return;
	pa = Kcalloc(km, mb128_t, n_hit[0] + n_hit[1]);
	for (r = n_pa = 0; r < 2; ++r) {
		for (i = 0; i < n_hit[r]; ++i) {
			mb128_t *p = &pa[n_pa++];
			mb_hit_t *h = &hit[r][i];
			h->pp = 0;
			p->x = l2b->ctg[h->tid].off + (h->rev? h->te : h->ts);
			p->y = (uint64_t)i << 2 | (uint64_t)h->rev << 1 | r;
		}
	}
	radix_sort_mb128x(pa, pa + n_pa);

	y[0] = y[1] = y[2] = y[3] = -1;
	for (i = 0; i < n_pa; ++i) {
		mb128_t *pi = &pa[i];
		mb_hit_t *hi = &hit[pi->y&1][pi->y>>2];
		for (r = 0; r < 2; ++r) {
			int which, dir = r << 1 | (pi->y>>1&1);
			if (pes[dir].failed) continue; // invalid orientation
			which = r << 1 | ((pi->y&1) ^ 1);
			if (y[which] < 0) continue; // no previous hit
			for (k = y[which]; k >= 0; --k) {
				mb128_t *pk = &pa[k], *q;
				mb_hit_t *hk = &hit[pk->y&1][pk->y>>2];
				int64_t dist;
				double ns, s;
				if ((pk->y&3) != which) continue;
				if (hi->tid != hk->tid) break;
				dist = pi->x - pk->x;
				if (dist > pes[dir].hi) break;
				if (dist < pes[dir].lo) continue;
				hk->pp = hi->pp = 1; // paired
				ns = (dist - pes[dir].avg) / pes[dir].std; // normalized score
				s = hk->p->dp_max + hi->p->dp_max + .721 * log(2. * erfc(fabs(ns) * MB_SQRT1_2)) * opt->a; // .721 = 1/log(4)
				if (s < 0.) s = 0.;
				if (n_pp == m_pp) Kgrow(km, mb128_t, pp, n_pp, m_pp);
				q = &pp[n_pp++];
				q->y = (pk->y&1) == 0? (uint64_t)k << 32 | i : (uint64_t)i << 32 | k; // upper bits: index to read[0]; lower: to read[1]
				q->x = (uint64_t)(s + .499) << 32 | ((hk->hash ^ hi->hash) & 0xffffffffULL);
			}
		}
		y[pi->y&3] = i;
	}

	ret->n_pp = n_pp;
	if (n_pp > 0) {
		uint64_t max = 0, max2 = 0;
		int32_t tmp = opt->a + opt->b > opt->q + opt->e? opt->a + opt->b : opt->q + opt->e;
		for (i = 0; i < n_pp; ++i) { // find max and max2
			mb128_t *q = &pp[i];
			if (q->x > max) max = q->x, ret->i[0] = q->y>>32, ret->i[1] = (uint32_t)q->y;
			else if (q->x > max2) max2 = q->x;
		}
		ret->score = max>>32, ret->sub_sc = max2>>32;
		for (i = 0; i < n_pp; ++i)
			if (pp[i].x>>32 <= ret->score - tmp)
				ret->n_sub++;
	}
	kfree(km, pp);
	kfree(km, pa);
}

static void mb_matesw_align(void *km, const mb_opt_t *opt, int32_t qlen, uint8_t *qseq, int32_t tlen, uint8_t *tseq, mb_hit_t *h, ksw_extz_t *ez)
{
	int8_t mat[25];
	int32_t max_sc = qlen < tlen? qlen : tlen;
	int32_t b_mm = (opt->b + opt->a - 1) / opt->a;
	int32_t b_ts = (opt->b_ts + opt->a - 1) / opt->a;
	int32_t b_ambi = (opt->b_ambi + opt->a - 1) / opt->a;
	int32_t gapo = (opt->q + opt->a - 1) / opt->a;
	int32_t gape = (opt->e + opt->a - 1) / opt->a;
	int32_t sz, xtra;
	void *qp;
	ksw_llrst_t rst;

	memset(h, 0, sizeof(*h));
	if (max_sc >= 32767) return;
	ksw_gen_ts_mat(5, mat, 1, b_mm, b_ts, b_ambi);
	sz = max_sc < 255 - b_mm? 1 : 2;
	qp = ksw_ll_qinit(km, sz, qlen, qseq, 5, mat);
	xtra = KSW_LL_SUBO | opt->min_len;
	if (sz == 1)
		rst = ksw_ll_u8_core(qp, tlen, tseq, gapo, gape, xtra);
	else
		rst = ksw_ll_i16_core(qp, tlen, tseq, gapo, gape, xtra);
	kfree(km, qp);
	if (rst.score * opt->a >= opt->min_dp_max) {
		int32_t te = rst.te + 1, qe = rst.qe + 1;
		mb_seq_rev(qe, qseq);
		mb_seq_rev(te, tseq);
		ksw_gen_ts_mat(5, mat, opt->a, opt->b, opt->b_ts, opt->b_ambi);
		ksw_extz2_sse(km, qe, qseq, te, tseq, 5, mat, opt->q, opt->e, opt->bw, opt->zdrop, opt->end_bonus, KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
		if (ez->n_cigar > 0 && ez->max >= opt->min_dp_max) {
			mb_append_cigar(h, ez->n_cigar, ez->cigar);
			h->qe = qe, h->te = te;
			h->ts = te - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
			h->qs = qe - (ez->reach_end? 0 : ez->max_q + 1);
			h->p->dp_max = ez->max;
			h->p->dp_score = ez->score;
			h->p->dp_max2 = (int32_t)((double)ez->max / rst.score * rst.score2 + .499);
			h->score = h->score0 = rst.score;
			h->subsc = rst.score2;
			h->cnt = 0, h->as = -1;
			h->parent = MB_PARENT_UNSET;
			kfree(km, ez->cigar);
		}
		mb_seq_rev(qe, qseq);
		mb_seq_rev(te, tseq);
	}
}

typedef struct {
	int32_t n, m;
	mb_hit_t *a;
} mb_hit_v;

static void mb_matesw_core(void *km, const mb_opt_t *opt, const l2b_t *l2b, const mb_pestat_t pes[4], const mb_hit_t *h0, int32_t len, uint8_t *seq[2], mb_hit_v *h1, ksw_extz_t *ez)
{
	int32_t r, skip[4];
	int64_t pos5;
	// find permitted orientation
	for (r = 0; r < 4; ++r) skip[r] = !!pes[r].failed;
	if (skip[0] + skip[1] + skip[2] + skip[3] == 4) return; // no need to perform SW
	// perform SW
	pos5 = h0->rev? h0->te : h0->ts;
	for (r = 0; r < 4; ++r) {
		int is_rev, is_large;
		int64_t ts, te;
		if (skip[r]) continue;
		is_rev = (r>>1 != (r&1)); // whether to reverse complement the mate
		is_large = !(r>>1); // whether the mate has larger coordinate
		ts = (is_large? pos5 + pes[r].lo : pos5 - pes[r].hi) - (!is_rev? 0 : len);
		te = (is_large? pos5 + pes[r].hi : pos5 - pes[r].lo) + (!is_rev? len : 0);
		if (ts < 0) ts = 0;
		if (te > l2b->ctg[h0->tid].len) te = l2b->ctg[h0->tid].len;
		if (te - ts > len) {
			uint8_t *ref;
			mb_hit_t ht;
			ref = Kmalloc(km, uint8_t, te - ts);
			l2b_getseq(l2b, h0->tid, ts, te, ref);
			mb_matesw_align(km, opt, len, seq[is_rev], te - ts, ref, &ht, ez);
			if (ht.p) { // a good hit found
				ht.tid = h0->tid;
				ht.ts += ts, ht.te += ts;
				ht.rev = is_rev;
				if (is_rev) {
					int32_t qt = ht.qs;
					ht.qs = len - ht.qe;
					ht.qe = len - qt;
				}
				if (h1->n == h1->m) kom_grow(mb_hit_t, h1->a, h1->n, h1->m);
				h1->a[h1->n++] = ht;
			}
			kfree(km, ref);
		}
	}
}

static int32_t mb_matesw(void *km, const mb_opt_t *opt, const l2b_t *l2b, int32_t n_hit[2], mb_hit_t *hit[2], const mb_pestat_t pes[4], int32_t qlen[2], char *const qseq[2])
{
	int32_t r, i, j, n_res, max[2] = {0, 0}, ori_sum = n_hit[0] + n_hit[1];
	uint64_t rng = 0, *a;
	mb_hit_v ha[2];
	uint8_t *qs[2][2];
	ksw_extz_t ez;

	// collect rescue candidates
	if (opt->max_rescue == 0) return 0;
	for (r = 0; r < 2; ++r) // find the max score
		for (i = 0; i < n_hit[r]; ++i)
			rng ^= hit[r][i].hash, max[r] = max[r] > hit[r][i].p->dp_max? max[r] : hit[r][i].p->dp_max;
	for (r = 0, n_res = 0; r < 2; ++r)
		for (i = 0; i < n_hit[r]; ++i)
			if (hit[r][i].pp == 0 && hit[r][i].p->dp_max >= max[r] - opt->pen_unpair)
				++n_res;
	if (n_res == 0) return 0; // nothing to rescue
	if (n_res > opt->max_rescue) n_res = opt->max_rescue;
	a = Kcalloc(km, uint64_t, n_res);
	for (r = j = 0; r < 2; ++r) { // candidates for rescue
		for (i = 0; i < n_hit[r]; ++i) {
			if (hit[r][i].pp == 0 && hit[r][i].p->dp_max >= max[r] - opt->pen_unpair) { // reservior sampling
				int32_t y;
				y = j++ < n_res? j - 1 : (int32_t)(j * kom_u64todbl(kom_splitmix64(&rng)));
				if (y < n_res) a[y] = (uint64_t)r << 32 | i;
			}
		}
		ha[r].n = ha[r].m = n_hit[r], ha[r].a = hit[r];
	}

	// prepare query sequences
	qs[0][0] = Kcalloc(km, uint8_t, (qlen[0] + qlen[1]) * 2);
	qs[0][1] = qs[0][0] + qlen[0];
	qs[1][0] = qs[0][1] + qlen[0];
	qs[1][1] = qs[1][0] + qlen[1];
	for (r = 0; r < 2; ++r) {
		for (i = 0; i < qlen[r]; ++i) {
			int32_t c = kom_nt4_table[(uint8_t)qseq[r][i]];
			qs[r][0][i] = c;
			qs[r][1][qlen[r] - 1 - r] = c < 4? 3 - c : 4;
		}
	}

	// do alignment
	memset(&ez, 0, sizeof(ez));
	ez.m_cigar = 8;
	ez.cigar = Kmalloc(km, uint32_t, ez.m_cigar);
	for (i = 0; i < n_res; ++i) {
		int32_t r = a[i]>>32&1, j = (int32_t)a[i];
		mb_matesw_core(km, opt, l2b, pes, &hit[r][j], qlen[!r], qs[!r], &ha[!r], &ez);
	}
	for (r = 0; r < 2; ++r)
		n_hit[r] = ha[r].n, hit[r] = ha[r].a;

	kfree(km, ez.cigar);
	kfree(km, qs[0][0]);
	kfree(km, a);
	return n_hit[0] + n_hit[1] - ori_sum;
}

void mb_pair(void *km, const mb_opt_t *opt, const l2b_t *l2b, int32_t n_hit[2], mb_hit_t *hit[2], const mb_pestat_t pes[4], int32_t qlen[2], char *const qseq[2])
{
	mb_pairaux_t paux;
	mb_pair_hits(km, opt, l2b, n_hit, hit, pes, &paux);
	if (opt->max_rescue > 0) {
		if (mb_matesw(km, opt, l2b, n_hit, hit, pes, qlen, qseq) > 0)
			mb_pair_hits(km, opt, l2b, n_hit, hit, pes, &paux); // pair again if new hits rescued
	}
}
