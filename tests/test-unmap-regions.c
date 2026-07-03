#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "mbpriv.h"

static mb_hit_t hit(int tid, int64_t st, int64_t en)
{
	mb_hit_t h;
	memset(&h, 0, sizeof(h));
	h.tid = tid;
	h.ts = st;
	h.te = en;
	return h;
}

static int32_t depth(const mb_unmap_regions_t *r, int tid, int64_t st, int64_t en)
{
	mb_hit_t h = hit(tid, st, en);
	return mb_unmap_hit_max_depth(r, &h);
}

int main(int argc, char **argv)
{
	l2b_t l2b;
	l2b_ctg_t ctg[2];
	mb_unmap_regions_t *r;
	mb_hit_t h[2];
	if (argc != 2) return 2;
	memset(&l2b, 0, sizeof(l2b));
	memset(ctg, 0, sizeof(ctg));
	ctg[0].name = "chr1";
	ctg[0].len = 1000;
	ctg[1].name = "2";
	ctg[1].len = 1000;
	ctg[1].off = 1000;
	l2b.n_ctg = 2;
	l2b.ctg = ctg;

	r = mb_unmap_regions_load(argv[1], &l2b);
	assert(r);
	assert(r->n_regions == 4);

	assert(depth(r, 0, 99, 100) == 0);
	assert(depth(r, 0, 100, 101) == 40);
	assert(depth(r, 0, 199, 200) == 40);
	assert(depth(r, 0, 200, 201) == 0);
	assert(depth(r, 1, 9, 10) == 7);
	assert(depth(r, 1, 20, 21) == 0);

	assert(depth(r, 0, 500, 501) == 9);
	assert(depth(r, 0, 301, 302) == 9);
	assert(depth(r, 0, 302, 303) == 9);
	assert(depth(r, 0, 300, 301) == 9);
	assert(depth(r, 0, 799, 800) == 9);
	assert(depth(r, 0, 299, 300) == 0);
	assert(depth(r, 0, 800, 801) == 0);

	h[0] = hit(0, 100, 105);
	h[1] = hit(0, 200, 205);
	mb_apply_unmap_regions(r, 2, h);
	assert(h[0].unmap == 1 && h[0].unmap_max_depth == 40);
	assert(h[1].unmap == 0 && h[1].unmap_max_depth == 0);

	mb_unmap_regions_destroy(r);
	puts("ok");
	return 0;
}
