#include <stdlib.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "kommon.h"
#include "ksort.h"

#define key_128x(a) ((a).x)
KRADIX_SORT_INIT(mb128x, mb128_t, key_128x, 8)

struct mb_tbuf_s {
	void *km;
};

mb_tbuf_t *mb_tbuf_init(void)
{
	mb_tbuf_t *b;
	b = kom_calloc(mb_tbuf_t, 1);
	if (!(kom_dbg_flag & MB_DBG_NO_KALLOC))
		b->km = km_init();
	return b;
}

void mb_tbuf_destroy(mb_tbuf_t *b)
{
	if (b->km) km_destroy(b->km);
	free(b);
}

static inline uint64_t hash64(uint64_t x)
{
	x ^= x >> 30;
	x *= 0xbf58476d1ce4e5b9ULL;
	x ^= x >> 27;
	x *= 0x94d049bb133111ebULL;
	x ^= x >> 31;
	return x;
}

static inline void mb_hit_set_coor(mb_hit_t *r, int32_t qlen, const l2b_t *l2b, const mb_anchor_t *a)
{ // NB: r->as and r->cnt MUST BE set correctly for this function to work
	int32_t k = r->as;
	const mb_anchor_t *ak = &a[k];
	const mb_anchor_t *ak_last = &a[k + r->cnt - 1];

	if (ak->tid < l2b->n_ctg) {
		r->tid = ak->tid, r->rev = 0;
		r->ts = ak->pos - l2b->ctg[r->tid].off;
		r->te = ak_last->pos + ak_last->len - l2b->ctg[r->tid].off;
	} else {
		r->tid = 2 * l2b->n_ctg - 1 - ak->tid, r->rev = 1;
		r->ts = l2b->tot_len * 2 - (ak_last->pos + ak_last->len) - l2b->ctg[r->tid].off;
		r->te = l2b->tot_len * 2 - ak->pos - l2b->ctg[r->tid].off;
	}
	r->qs = ak->qs;
	r->qe = ak_last->qs + ak_last->len;
}

mb_hit_t *mb_gen_hit(void *km, uint32_t hash, int qlen, const mb_idx_t *idx, int n_u, uint64_t *u, mb_anchor_t *a)
{ // convert chains to hits
	mb128_t *z, tmp;
	mb_hit_t *r;
	int i, k;

	if (n_u <= 0) return 0;

	// sort by score
	z = Kmalloc(km, mb128_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		uint32_t h;
		h = (uint32_t)hash64((hash64(a[k].pos) + hash64(a[k].qs)) ^ hash);
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
		ri->score = z[i].x >> 32;
		ri->cnt = (int32_t)z[i].y;
		ri->as = z[i].y >> 32;
		mb_hit_set_coor(ri, qlen, idx->l2b, a);
	}
	kfree(km, z);
	return r;
}

mb_hit_t *mb_map(const mb_idx_t *idx, int64_t qlen, const char *seq0, int32_t *n_hit, mb_tbuf_t *b, const mb_mopt_t *opt, const char *qname)
{
	int64_t i;
	void *km = b->km;
	uint8_t *seq;

	seq = Kcalloc(km, uint8_t, qlen);
	for (i = 0; i < qlen; ++i)
		seq[i] = kom_nt4_table[(uint8_t)seq0[i]];
	kfree(km, seq);
	return 0;
}
