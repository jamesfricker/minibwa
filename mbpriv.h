#ifndef MBPRIV_H
#define MBPRIV_H

#include "minibwa.h"
#include "l2bit.h"
#include "bwt.h"

#define MB_DBG_NO_KALLOC    (0x1LL)

struct mb_idx_s {
	l2b_t *l2b;
	mb_bwt_t *bwt;
};

typedef struct {
	int32_t tid, len;
	int32_t qs;
	uint16_t qocc, tocc:15, flt:1;
	uint64_t pos;
} mb_anchor_t;

typedef struct { int64_t n, m; mb_anchor_t *a; } mb_anchor_v;

typedef struct { uint64_t x, y; } mb128_t;
void radix_sort_mb128x(mb128_t *beg, mb128_t *end);

#ifdef __cplusplus
extern "C" {
#endif

// defined in bwtgen.c
void mb_bwtgen(const char *fn_pac, const char *fn_bwt, int block_size);

// defined in seed.c
void mb_seed_intv(void *km, const mb_bwt_t *bwt, int32_t len, const uint8_t *seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v);
void mb_anchor(void *km, const mb_idx_t *idx, const mb_sai_v *u, int32_t max_occ, mb_anchor_v *v);
void mb_anchor_sort(int64_t n, mb_anchor_t *a);

// defined in lchain.c
mb_anchor_t *mb_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
						  int64_t n, mb_anchor_t *a, int *n_u_, uint64_t **_u, void *km);
mb_anchor_t *mb_lchain_rmq(int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
						   int64_t n, mb_anchor_t *a, int *n_u_, uint64_t **_u, void *km);

// defined in map-algo.c
mb_hit_t *mb_gen_hit(void *km, uint32_t hash, int qlen, const mb_idx_t *idx, int n_u, uint64_t *u, mb_anchor_t *a);

#ifdef __cplusplus
}
#endif

#endif
