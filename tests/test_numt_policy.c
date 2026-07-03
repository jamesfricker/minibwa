/* End-to-end style unit test for the --numt (MB_F_NUMT) chrM-vs-nuclear policy.
 *
 * It exercises the public policy entry points mb_apply_numt_primary() and
 * mb_apply_numt_mapq() together with the SAM formatter to demonstrate, in one
 * transcript, the behaviours described in the README/manpage:
 *
 *   1. exact score ties between chrM and nuclear hits deterministically prefer
 *      the nuclear placement (chrM hit is demoted to secondary);
 *   2. near-tie chrM/nuclear competitors are MAPQ-capped (to opt.numt_mapq_cap)
 *      and annotated with ng:Z:chrM-nuclear;
 *   3. a read unique to chrM keeps its high MAPQ and is NOT tagged;
 *   4. with the policy disabled (no MB_F_NUMT) nothing changes.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mbpriv.h"

static mb_extra_t *make_extra(int dp_max)
{
	mb_extra_t *p = (mb_extra_t*)calloc(1, sizeof(*p) + sizeof(uint32_t));
	p->cap = 1;
	p->dp_score = dp_max;
	p->dp_max0 = dp_max;
	p->dp_max = dp_max;
	p->dp_max2 = 0;
	p->n_cigar = 1;
	p->cigar[0] = 150 << 4 | MB_CIGAR_MATCH;
	return p;
}

static mb_hit_t make_hit(int tid, int id, int mapq, int dp_max)
{
	mb_hit_t h;
	memset(&h, 0, sizeof(h));
	h.tid = tid;
	h.ts = 100;
	h.te = 250;
	h.id = id;
	h.parent = id;          /* each hit starts as its own primary */
	h.score = dp_max;
	h.qs = 0;
	h.qe = 150;
	h.mlen = 150;
	h.blen = 150;
	h.mapq = mapq;
	h.p = make_extra(dp_max);
	return h;
}

static void free_hits(mb_hit_t *h, int n)
{
	int i;
	for (i = 0; i < n; ++i) free(h[i].p);
}

static int primary_tid(int n, const mb_hit_t *h)
{
	int i;
	for (i = 0; i < n; ++i)
		if (h[i].parent == h[i].id) return (int)h[i].tid;
	return -1;
}

static const mb_hit_t *find_tid(int n, const mb_hit_t *h, int tid)
{
	int i;
	for (i = 0; i < n; ++i)
		if ((int)h[i].tid == tid) return &h[i];
	return NULL;
}

int main(void)
{
	l2b_t l2b;
	l2b_ctg_t ctg[2];
	mb_opt_t opt;

	memset(&l2b, 0, sizeof(l2b));
	memset(ctg, 0, sizeof(ctg));
	ctg[0].name = "chr1"; ctg[0].len = 100000;  /* nuclear */
	ctg[1].name = "chrM"; ctg[1].len = 16569;   /* mitochondrial */
	l2b.n_ctg = 2;
	l2b.ctg = ctg;

	mb_opt_init(&opt);
	assert(opt.numt_mapq_cap == 10);           /* documented default cap */

	/* ---- Case A: exact score tie, chrM(primary) vs nuclear ---- */
	{
		mb_hit_t h[2];
		int promoted;
		h[0] = make_hit(1, 0, 60, 300);   /* chrM, currently primary, high mapq */
		h[1] = make_hit(0, 1, 0, 300);    /* nuclear, same DP score (exact tie) */
		h[1].parent = h[0].id;            /* nuclear is a secondary of the chrM primary */

		opt.flag |= MB_F_NUMT;
		promoted = mb_apply_numt_primary(&opt, &l2b, 2, h);
		assert(promoted == 1);
		assert(primary_tid(2, h) == 0);   /* nuclear (chr1) is now the primary */
		printf("[A] exact tie -> primary tid=%d (chr1=nuclear), chrM demoted\n", primary_tid(2, h));

		mb_apply_numt_mapq(&opt, &l2b, 2, h);
		assert(h[0].numt_ambig == 1 && h[1].numt_ambig == 1);
		free_hits(h, 2);
	}

	/* ---- Case B: near tie -> MAPQ capped to 10 and ng:Z tagged ---- */
	{
		mb_hit_t h[2];
		int mt_dp = 300, nuc_dp;
		nuc_dp = mt_dp - opt.numt_score_diff;  /* exactly at the ambiguity threshold */
		h[0] = make_hit(1, 0, 60, mt_dp);      /* chrM primary, naturally MAPQ 60 */
		h[1] = make_hit(0, 1, 60, nuc_dp);     /* nuclear, near-tie score */
		h[1].parent = h[0].id;

		mb_apply_numt_mapq(&opt, &l2b, 2, h);
		assert(h[0].numt_ambig == 1 && h[1].numt_ambig == 1);
		assert(h[0].mapq <= opt.numt_mapq_cap);  /* chrM primary capped */
		printf("[B] near tie -> chrM MAPQ capped 60 -> %d, both ng-tagged\n", h[0].mapq);
		free_hits(h, 2);
	}

	/* ---- Case C: read unique to chrM keeps high MAPQ, not tagged ---- */
	{
		mb_hit_t h[1];
		h[0] = make_hit(1, 0, 60, 300);   /* chrM only */
		mb_apply_numt_primary(&opt, &l2b, 1, h);
		mb_apply_numt_mapq(&opt, &l2b, 1, h);
		assert(h[0].numt_ambig == 0);
		assert(h[0].mapq == 60);
		printf("[C] unique chrM read keeps MAPQ=%d, ng untagged\n", h[0].mapq);
		free_hits(h, 1);
	}

	/* ---- Case D: policy OFF leaves chrM primary + high MAPQ untouched ---- */
	{
		mb_hit_t h[2];
		int promoted;
		mb_opt_t off;
		mb_opt_init(&off);                /* MB_F_NUMT not set */
		h[0] = make_hit(1, 0, 60, 300);   /* chrM primary */
		h[1] = make_hit(0, 1, 0, 300);    /* nuclear, exact tie */
		h[1].parent = h[0].id;
		promoted = mb_apply_numt_primary(&off, &l2b, 2, h);
		mb_apply_numt_mapq(&off, &l2b, 2, h);
		assert(promoted == 0);
		assert(primary_tid(2, h) == 1);   /* still chrM */
		assert(h[0].numt_ambig == 0 && h[1].numt_ambig == 0);
		assert(h[0].mapq == 60);
		printf("[D] policy OFF -> primary tid=%d (chrM), MAPQ=%d, untagged\n",
			primary_tid(2, h), h[0].mapq);
		free_hits(h, 2);
	}

	/* ---- Case E: formatter emits ng:Z:chrM-nuclear only when ambiguous ---- */
	{
		mb_bseq1_t read;
		mb_hit_t h[2], *hitp[1];
		int32_t n = 2, i;
		char seq[151];
		kstring_t out = {0,0,0};
		const mb_hit_t *chrm;

		for (i = 0; i < 150; ++i) seq[i] = "ACGT"[i & 3];
		seq[150] = 0;
		memset(&read, 0, sizeof(read));
		read.name = "read_numt";
		read.seq = seq;
		read.l_seq = 150;   /* matches the 150M CIGAR built by make_extra() */

		h[0] = make_hit(0, 0, 10, 300);   /* nuclear primary after promotion */
		h[1] = make_hit(1, 1, 0, 300);    /* chrM secondary */
		h[1].parent = h[0].id;
		h[0].sam_pri = 1;

		opt.flag = MB_F_NUMT;
		opt.xa_max = 0;
		mb_apply_numt_mapq(&opt, &l2b, 2, h);
		chrm = find_tid(2, h, 1);
		assert(chrm && chrm->numt_ambig == 1);

		hitp[0] = h;
		mb_format(0, &out, &l2b, &read, 1, &n, hitp, 0, &opt, 0, 0);
		assert(strstr(out.s, "ng:Z:chrM-nuclear") != NULL);
		printf("[E] SAM line carries ng:Z:chrM-nuclear tag\n");
		free(out.s);
		free_hits(h, 2);
	}

	printf("ALL NUMT POLICY CHECKS PASSED\n");
	return 0;
}
