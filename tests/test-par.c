#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mbpriv.h"

static mb_extra_t *make_extra(int32_t dp)
{
	mb_extra_t *p = (mb_extra_t*)calloc(1, sizeof(mb_extra_t));
	assert(p);
	p->dp_score = dp;
	p->dp_max0 = dp;
	p->dp_max = dp;
	return p;
}

int main(void)
{
	l2b_ctg_t ctg[2];
	l2b_t l2b;
	mb_hit_t h[2];
	mb_bseq1_t read;
	kstring_t paf = {0, 0, 0};

	memset(&l2b, 0, sizeof(l2b));
	memset(ctg, 0, sizeof(ctg));
	ctg[0].name = "chrX";
	ctg[0].len = 156040895;
	ctg[1].name = "chrY";
	ctg[1].len = 57227415;
	l2b.n_ctg = 2;
	l2b.ctg = ctg;
	mb_par_init(&l2b);
	assert(l2b.par_asm == 38);

	memset(h, 0, sizeof(h));
	h[0].tid = 1; /* chrY sorts first in this fixture to exercise primary promotion. */
	h[0].ts = 11000;
	h[0].te = 11100;
	h[0].qs = 0;
	h[0].qe = 100;
	h[0].score = h[0].score0 = 100;
	h[0].mlen = h[0].blen = 100;
	h[0].seed_ratio = 255;
	h[0].p = make_extra(200);

	h[1] = h[0];
	h[1].tid = 0;
	h[1].p = make_extra(200);

	mb_mark_par_hits(&l2b, 2, h);
	assert(h[0].par == 1);
	assert(h[1].par == 1);
	assert(mb_par_equiv(&l2b, &h[0], &h[1]));

	mb_set_parent(0, &l2b, 0.5f, 0x7fffffff, 2, h, 10, 0);
	assert(h[0].parent == 0);
	assert(h[1].parent == 0);
	assert(h[0].subsc == 0);
	assert(h[0].p->dp_max2 == 0);

	mb_par_resolve(&l2b, 2, h, 10);
	assert(h[1].parent == 1);
	assert(h[0].parent == 1);
	mb_set_sam_pri(2, h, 0);
	assert(h[1].sam_pri == 1);
	assert(h[0].sam_pri == 0);

	mb_set_mapq(0, 100, 2, h, 25, 2, 1, 325);
	assert(h[1].mapq >= 10);
	assert(h[0].mapq == 0);

	memset(&read, 0, sizeof(read));
	read.name = "par-read";
	read.l_seq = 100;
	mb_fmt_paf(&paf, &l2b, &read, &h[1], 0, 1, 0);
	assert(strstr(paf.s, "\tpa:Z:PAR1"));

	printf("PAR regression ok: primary=chrX mapq=%d tag=%s\n", h[1].mapq, mb_par_name(h[1].par));
	free(paf.s);
	free(h[0].p);
	free(h[1].p);
	return 0;
}
