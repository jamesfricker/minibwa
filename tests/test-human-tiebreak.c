#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "mbpriv.h"

static mb_extra_t *make_extra(int32_t dp)
{
	mb_extra_t *p = (mb_extra_t*)calloc(1, sizeof(*p));
	assert(p);
	p->dp_max = dp;
	p->dp_max0 = dp;
	p->dp_score = dp;
	return p;
}

static mb_hit_t make_hit(int32_t tid, int32_t score)
{
	mb_hit_t h;
	memset(&h, 0, sizeof(h));
	h.tid = tid;
	h.ts = 100;
	h.te = 110;
	h.qs = 0;
	h.qe = 10;
	h.score = h.score0 = score;
	h.mlen = h.blen = 10;
	h.parent = MB_PARENT_UNSET;
	h.p = make_extra(score);
	return h;
}

static void free_hits(mb_hit_t *h, int32_t n)
{
	int32_t i;
	for (i = 0; i < n; ++i)
		free(h[i].p);
}

int main(void)
{
	l2b_t l2b;
	l2b_ctg_t ctg[5];
	mb_hit_t h[4];

	memset(&l2b, 0, sizeof(l2b));
	memset(ctg, 0, sizeof(ctg));
	ctg[0].name = "HLA-A*01:01:01:01";
	ctg[1].name = "chr1_KI270762v1_alt";
	ctg[2].name = "hs37d5";
	ctg[3].name = "chr1";
	ctg[4].name = "NC_001422.1_phix";
	l2b.n_ctg = 5;
	l2b.ctg = ctg;

	h[0] = make_hit(0, 100);
	h[1] = make_hit(1, 100);
	h[2] = make_hit(2, 100);
	h[3] = make_hit(3, 100);
	mb_human_sort_ties(0, &l2b, 4, h);
	assert(h[0].tid == 3);
	assert(h[1].tid == 1);
	assert(h[2].tid == 0);
	assert(h[3].tid == 2);

	mb_mark_par_hits(&l2b, 4, h);
	mb_set_parent(0, &l2b, 0.5f, 0x7fffffff, 4, h, 10, 0);
	mb_set_sam_pri(4, h, 0);
	assert(h[0].tid == 3);
	assert(h[0].parent == 0);
	assert(h[0].sam_pri == 1);
	assert(h[1].parent == 0);
	assert(h[2].parent == 0);
	assert(h[3].parent == 0);
	free_hits(h, 4);

	h[0] = make_hit(0, 101);
	h[1] = make_hit(3, 100);
	mb_human_sort_ties(0, &l2b, 2, h);
	assert(h[0].tid == 0);
	assert(h[1].tid == 3);
	free_hits(h, 2);

	h[0] = make_hit(4, 100);
	h[1] = make_hit(2, 100);
	mb_human_sort_ties(0, &l2b, 2, h);
	assert(h[0].tid == 2);
	assert(h[1].tid == 4);
	free_hits(h, 2);

	h[0] = make_hit(0, 100);
	h[1] = make_hit(3, 100);
	h[2] = make_hit(1, 100);
	mb_human_sort_ties(0, &l2b, 3, h);
	assert(h[0].tid == 3);
	assert(h[1].tid == 1);
	assert(h[2].tid == 0);
	free_hits(h, 3);

	return 0;
}
