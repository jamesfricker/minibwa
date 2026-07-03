#include <assert.h>
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
	p->n_cigar = 1;
	p->cigar[0] = 10 << 4 | MB_CIGAR_MATCH;
	return p;
}

static mb_hit_t make_hit(int tid, int id, int parent, int sam_pri, int dp_max)
{
	mb_hit_t h;
	memset(&h, 0, sizeof(h));
	h.tid = tid;
	h.ts = 100 + tid * 100;
	h.te = h.ts + 10;
	h.id = id;
	h.parent = parent;
	h.sam_pri = sam_pri;
	h.qs = 0;
	h.qe = 10;
	h.mlen = 10;
	h.blen = 10;
	h.mapq = sam_pri? 20 : 0;
	h.p = make_extra(dp_max);
	return h;
}

static void free_hits(mb_hit_t *h, int n)
{
	int i;
	for (i = 0; i < n; ++i)
		free(h[i].p);
}

static int has_tag_value(const char *sam, const char *needle)
{
	return strstr(sam, needle) != 0;
}

static int primary_tid(int n, const mb_hit_t *h)
{
	int i;
	for (i = 0; i < n; ++i)
		if (h[i].parent == h[i].id) return (int)h[i].tid;
	return -1;
}

int main(void)
{
	l2b_t l2b;
	l2b_ctg_t ctg[4];
	mb_bseq1_t read;
	mb_opt_t opt;
	int32_t n_hit = 4;
	mb_hit_t hits[4], *hitp[1];
	kstring_t out = {0,0,0};

	memset(&l2b, 0, sizeof(l2b));
	memset(ctg, 0, sizeof(ctg));
	ctg[0].name = "chr6";
	ctg[0].len = 1000;
	ctg[1].name = "chr6_GL000250v2_alt";
	ctg[1].len = 1000;
	ctg[2].name = "HLA-A*01:01";
	ctg[2].len = 1000;
	ctg[3].name = "chr2";
	ctg[3].len = 1000;
	l2b.n_ctg = 4;
	l2b.ctg = ctg;

	memset(&read, 0, sizeof(read));
	read.name = "r1";
	read.seq = "ACGTACGTAA";
	read.l_seq = 10;

	hits[0] = make_hit(0, 0, 0, 1, 20);
	hits[1] = make_hit(1, 1, 1, 0, 19);
	hits[2] = make_hit(2, 2, 2, 0, 18);
	hits[3] = make_hit(3, 3, 3, 0, 18);
	hitp[0] = hits;

	mb_opt_init(&opt);
	opt.flag = MB_F_HUMAN_ALT;
	opt.xa_max = 5;
	opt.xa_ratio = 0.8f;
	mb_format(0, &out, &l2b, &read, 1, &n_hit, hitp, 0, &opt, 0, 0);
	assert(has_tag_value(out.s, "\tn2:i:2"));
	assert(has_tag_value(out.s, "\tXA:Z:chr6_GL000250v2_alt,+201,10M,0;HLA-A*01:01,+301,10M,0;"));
	assert(!has_tag_value(out.s, "SA:Z:chr6_GL000250v2_alt"));
	assert(!has_tag_value(out.s, "SA:Z:HLA-A*01:01"));
	assert(has_tag_value(out.s, "SA:Z:chr2,401,+,10M,0,0;"));

	out.l = 0;
	opt.xa_max = 1;
	mb_format(0, &out, &l2b, &read, 1, &n_hit, hitp, 0, &opt, 0, 0);
	assert(has_tag_value(out.s, "\tn2:i:2"));
	assert(!has_tag_value(out.s, "\tXA:Z:"));

	out.l = 0;
	opt.xa_max = 0;
	mb_format(0, &out, &l2b, &read, 1, &n_hit, hitp, 0, &opt, 0, 0);
	assert(!has_tag_value(out.s, "\tn2:i:"));
	assert(!has_tag_value(out.s, "\tXA:Z:"));
	assert(has_tag_value(out.s, "SA:Z:chr6_GL000250v2_alt"));
	assert(has_tag_value(out.s, "HLA-A*01:01"));

	{
		mb_hit_t h[2];
		int promoted;
		h[0] = make_hit(0, 0, 1, 0, 30); /* chr6 HLA-region placement as secondary */
		h[1] = make_hit(2, 1, 1, 1, 30); /* HLA allele contig currently primary */
		opt.hla_policy = MB_HLA_POLICY_MAIN_CONTIG;
		promoted = mb_apply_hla_primary(&opt, &l2b, 2, h);
		assert(promoted == 1);
		assert(primary_tid(2, h) == 0);
		free_hits(h, 2);
	}

	{
		mb_hit_t h[2];
		int promoted;
		h[0] = make_hit(0, 0, 1, 0, 30);
		h[1] = make_hit(2, 1, 1, 1, 30);
		opt.hla_policy = MB_HLA_POLICY_ALLELE_CONTIG;
		promoted = mb_apply_hla_primary(&opt, &l2b, 2, h);
		assert(promoted == 0);
		assert(primary_tid(2, h) == 2);
		free_hits(h, 2);
	}

	{
		mb_hit_t h[2];
		int promoted;
		h[0] = make_hit(0, 0, 1, 0, 29);
		h[1] = make_hit(2, 1, 1, 1, 30);
		opt.hla_policy = MB_HLA_POLICY_MAIN_CONTIG;
		promoted = mb_apply_hla_primary(&opt, &l2b, 2, h);
		assert(promoted == 0);
		assert(primary_tid(2, h) == 2);
		free_hits(h, 2);
	}

	free(out.s);
	free_hits(hits, 4);
	return 0;
}
