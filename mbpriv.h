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
	int32_t tid2; // stranded target sequence ID, ranged from 0 to 2*l2b_t::n_ctg - 1
	int32_t len; // length of the anchor
	int32_t qpos; // the query coordinate of the last base in the anchor; the start base is qpos+1-len
	uint32_t dummy; // this field can be used for bit flags or small information
	uint64_t tpos; // the target coordinate of the last base in the anchor
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
mb_anchor_t *mb_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_sc, float chn_pen_gap, float chn_pen_skip,
						  int64_t n, mb_anchor_t *a, int *n_u_, uint64_t **_u, void *km);
mb_anchor_t *mb_lchain_rmq(int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, int min_sc, float chn_pen_gap, float chn_pen_skip,
						   int64_t n, mb_anchor_t *a, int *n_u_, uint64_t **_u, void *km);

// defined in map-algo.c
mb_hit_t *mb_gen_hit(void *km, uint32_t hash, int qlen, const mb_idx_t *idx, int n_u, uint64_t *u, mb_anchor_t *a);
void mb_set_parent(void *km, float mask_level, int mask_len, int n, mb_hit_t *r, int sub_diff, int hard_mask_level);
void mb_hit_sort(void *km, int *n_regs, mb_hit_t *r);
void mb_sync_hits(void *km, int n_regs, mb_hit_t *regs);
void mb_select_sub(void *km, float pri_ratio, int min_diff, int best_n, int *n_, mb_hit_t *r);
void mb_filter_hits(const mb_mopt_t *opt, int qlen, int *n_regs, mb_hit_t *regs);
int mb_squeeze_a(void *km, int n_regs, mb_hit_t *regs, mb_anchor_t *a);
void mb_split_hit(mb_hit_t *r, mb_hit_t *r2, int n, int qlen, mb_anchor_t *a, const l2b_t *l2b);

#ifdef __cplusplus
}
#endif

#endif
