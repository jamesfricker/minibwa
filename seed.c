#include <string.h>
#include <stdio.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "ksort.h"

#define key_sai0(a) ((a).x[0])
KRADIX_SORT_INIT(mb_sai0, mb_sai_t, key_sai0, 8)

#define key_sais(a) ((a).size)
KRADIX_SORT_INIT(mb_sais, mb_sai_t, key_sais, 8)

static void mb_bwt_extend_eq(const mb_bwt_t *bwt, int32_t len, const uint8_t *q, mb_sai_t *p)
{
	uint32_t st = p->info>>32, en = (uint32_t)p->info;
	int32_t j, c;
	mb_sai_t ok[4];
	for (j = en; j < len; ++j) { // extend forwardly
		if (q[j] > 3) break;
		c = 3 - q[j];
		mb_bwt_extend(bwt, p, ok, 0);
		if (ok[c].size < p->size) break;
		*p = ok[c];
	}
	en = j;
	for (j = st - 1; j >= 0; --j) { // extend backwardly
		if (q[j] > 3) break;
		c = q[j];
		mb_bwt_extend(bwt, p, ok, 1);
		if (ok[c].size < p->size) break;
		*p = ok[c];
	}
	st = j + 1;
	p->info = (uint64_t)st<<32 | en;
}

static int32_t mb_bwt_seed_greedy(const mb_bwt_t *bwt, int32_t len, const uint8_t *q, int32_t x, int32_t min_len, int32_t max_intv, mb_sai_t *p)
{ // inspired by PMID:21209072 but more crude
	int32_t i, c;
	mb_sai_t ik, ok[4];
	memset(p, 0, sizeof(*p));
	if (q[x] > 3) return x + 1;
	mb_bwt_set_intv(bwt, q[x], &ik); // the initial interval of a single base
	for (i = x + 1; i < len; ++i) { // forward search
		if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			mb_bwt_extend(bwt, &ik, ok, 0);
			if (ok[c].size < max_intv && i - x >= min_len) {
				*p = ok[c];
				p->info = (uint64_t)x<<32 | (i + 1);
				return i + 1;
			}
			ik = ok[c];
		} else return i + 1;
	}
	return len;
}

void mb_seed_intv(void *km, const mb_bwt_t *bwt, int32_t len, const uint8_t *seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v)
{
	int64_t x = 0, i, n_a0;
	mb_sai_t p;

	v->n = 0;
	do {
		x = mb_bwt_smem(bwt, len, seq, x, min_len, 1, 1, &p);
		if (p.size > 0) {
			Kgrow(km, mb_sai_t, v->a, v->n, v->m);
			v->a[v->n++] = p;
		}
	} while (x < len);

	n_a0 = v->n;
	for (i = 0; i < n_a0; ++i) {
		uint32_t st = v->a[i].info>>32, en = (uint32_t)v->a[i].info;
		if (en - st < min_len * 2 || v->a[i].size > max_sub_occ)
			continue;
		x = st;
		do {
			x = mb_bwt_seed_greedy(bwt, en, seq, x, min_len, max_sub_occ * 2, &p);
			if (p.size > v->a[i].size) {
				int32_t to_add = 1;
				if (v->n > 0 && p.size == v->a[v->n-1].size && (uint32_t)p.info <= (uint32_t)v->a[v->n-1].info)
					to_add = 0;
				if (to_add) {
					mb_bwt_extend_eq(bwt, len, seq, &p);
					Kgrow(km, mb_sai_t, v->a, v->n, v->m);
					v->a[v->n++] = p;
				}
			}
		} while (x < en);
	}
}

void mb_anchor(void *km, const mb_idx_t *idx, const mb_sai_v *u, int32_t max_occ, mb_anchor_v *v)
{
	int64_t i, i0, j, k;
	uint64_t *a;

	v->n = 0;
	if (u->n == 0) return; // no anchors

	// sort by ::x[0] and then by ::size
	if (u->n > 1) radix_sort_mb_sai0(u->a, u->a + u->n);
	for (i = 1, i0 = 0; i <= u->n; ++i) {
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0]) {
			if (i - i0 > 1)
				radix_sort_mb_sais(&u->a[i0], &u->a[i]);
			i0 = i;
		}
	}

	for (i = 0, k = 0; i < u->n; ++i) // calculate the size of v->a
		k += u->a[i].size < max_occ? u->a[i].size : max_occ;
	Kgrow(km, mb_anchor_t, v->a, k - 1, v->m); // preallocate

	a = Kmalloc(km, uint64_t, max_occ * 2);
	for (i = 1, i0 = 0; i <= u->n; ++i) { // a bit overkilling for short reads, but may be beneficial for long centromeric reads
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0] || u->a[i].size != u->a[i0].size) {
			const mb_sai_t *p = &u->a[i0];
			int32_t n = 0, step = p->size < max_occ? 1 : p->size / max_occ;
			for (j = 0; j < p->size; j += step)
				a[n++] = p->x[0] + j;
			mb_bwt_sa_batch(km, idx->bwt, n, a);
			for (k = 0; k < n; ++k) {
				for (j = i0; j < i; ++j) {
					int32_t qs = u->a[j].info>>32, qe = (int32_t)u->a[j].info;
					int32_t rev, len = qe - qs;
					int64_t tid, cst;
					mb_anchor_t *q;
					tid = l2b_intv2cid(idx->l2b, a[k], a[k] + len, &cst, &rev);
					if (tid < 0) continue;
					Kgrow(km, mb_anchor_t, v->a, v->n, v->m);
					q = &v->a[v->n++];
					q->tid = rev? idx->l2b->n_ctg * 2 - 1 - q->tid : q->tid;
					q->len = len;
					q->qs = qs;
					q->pos = a[k];
					q->qocc = i - i0 < UINT16_MAX? i - i0 : UINT16_MAX;
					q->tocc = p->size < UINT16_MAX? p->size : UINT16_MAX;
				}
			}
		}
	}
	kfree(km, a);
}
