#ifndef L2BIT_H
#define L2BIT_H

#include <stdint.h>

#define L2B_MAGIC "L2B\1"

typedef struct {
	char *name, *comm;
	uint64_t len, off;
} l2b_ctg_t;

typedef struct {
	uint64_t st, en;
} l2b_intv_t;

typedef struct {
	uint64_t tot_len;
	uint64_t n_ctg, m_ctg;
	l2b_ctg_t *ctg;
	uint64_t n_pac, m_pac;
	uint64_t n_ambi, m_ambi;
	uint64_t n_mask, m_mask;
	l2b_intv_t *ambi, *mask;
	uint64_t *pac;
	char *cat_name, *cat_comm;
} l2b_t;

#ifdef __cplusplus
extern "C" {
#endif

l2b_t *l2b_import(const char *fn, uint64_t seed);
l2b_t *l2b_load(const char *fn);
int l2b_save(const char *fn, const l2b_t *l2b);
void l2b_destroy(l2b_t *l2b);

#ifdef __cplusplus
}
#endif

#endif
