#include <stdlib.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "kommon.h"
#include "ksort.h"

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

#define key_128x(a) ((a).x)
KRADIX_SORT_INIT(mb128x, mb128_t, key_128x, 8)

static inline uint64_t hash64(uint64_t key)
{
	key = (~key + (key << 21));
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8));
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4));
	key = key ^ key >> 28;
	key = (key + (key << 31));
	return key;
}

static inline void mb_hit_set_coor(mb_hit_t *r, int32_t qlen, const mb_idx_t *idx, const mb_anchor_t *a)
{ // NB: r->as and r->cnt MUST BE set correctly for this function to work
	int32_t k = r->as;
	const mb_anchor_t *ak = &a[k];
	const mb_anchor_t *ak_last = &a[k + r->cnt - 1];
	int64_t cst, cid;
	int rev;

	// Use l2b_intv2cid to convert global coordinates to contig-relative coordinates
	cid = l2b_intv2cid(idx->l2b, ak->pos, ak->pos + ak->len, &cst, &rev);
	r->tid = cid;
	r->rev = rev;
	r->ts = cst;
	
	// For the end coordinate, we need the relative coordinate of the last anchor
	l2b_intv2cid(idx->l2b, ak_last->pos, ak_last->pos + ak_last->len, &cst, &rev);
	r->te = cst + ak_last->len; // end is one-past-the-last base on the forward strand

	// query coordinates
	r->qs = ak->qs + 1 > ak->len ? ak->qs + 1 - ak->len : 0;
	r->qe = ak_last->qs + 1;
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
		mb_hit_set_coor(ri, qlen, idx, a);
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
