#include <stdio.h>
#include <stdlib.h>
#include "mbpriv.h"
#include "kommon.h"

static char mb_rg_id[256];

/**************
 * PAF output *
 **************/

void mb_fmt_paf(kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, const mb_hit_t *p, uint64_t opt_flag, int n_seg, int seg_idx)
{
	kom_sprintf_lite(s, "%s", t->name);
	if (n_seg > 1 && seg_idx >= 0)
		kom_sprintf_lite(s, "/%d", seg_idx + 1);
	kom_sprintf_lite(s, "\t%ld", (long)t->l_seq);
	if (p == 0) { // for unmapped reads
		kom_sprintf_lite(s, "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\n");
		return;
	}
	kom_sprintf_lite(s, "\t%d\t%d\t%c\t%s\t%ld\t%ld\t%ld\t%d\t%d\t%d\ttp:A:%c\ts1:i:%d\tcm:i:%d",
		p->qs, p->qe, p->rev? '-' : '+', l2b->ctg[p->tid].name, (long)l2b->ctg[p->tid].len, (long)p->ts, (long)p->te,
		p->mlen, p->blen, p->mapq, p->parent == p->id? 'P' : 'S', p->score, p->cnt);
	if (p->parent == p->id) kom_sprintf_lite(s, "\ts2:i:%d", p->subsc);
	if (p->p) {
		int32_t nm = p->blen - p->mlen + p->p->n_ambi;
		kom_sprintf_lite(s, "\tNM:i:%d\tAS:i:%d\tms:i:%d\tm2:i:%d", nm, p->p->dp_score, p->p->dp_max, p->p->dp_max2);
		if (p->p->n_cigar > 0) {
			int32_t i;
			kom_sprintf_lite(s, "\tcg:Z:");
			for (i = 0; i < p->p->n_cigar; ++i)
				kom_sprintf_lite(s, "%d%c", p->p->cigar[i]>>4, MB_CIGAR_STR[p->p->cigar[i]&0xf]);
		}
		if (opt_flag & (MB_F_WRITE_DS|MB_F_WRITE_CS|MB_F_WRITE_MD))
			kom_sprintf_lite(s, "\t%s", (char*)&p->p->cigar[p->p->n_cigar]);
	}
	if ((opt_flag & MB_F_COPY_COMMENT) && t->comment)
		kom_sprintf_lite(s, "\t%s", t->comment);
	kom_sprintf_lite(s, "\n");
}

/**************
 * SAM header *
 **************/

static char *mb_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}

static int sam_write_rg_line(kstring_t *str, const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(mb_rg_id, 0, 256);
	if (s == 0) return 0;
	if (strstr(s, "@RG") != s) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] the read group line is not started with @RG\n");
		goto err_set_rg;
	}
	if (strstr(s, "\t") != NULL) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n");
		goto err_set_rg;
	}
	rg_line = kom_strdup(s);
	mb_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] no ID within the read group line\n");
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] @RG:ID is longer than 255 characters\n");
		goto err_set_rg;
	}
	for (q = p, r = mb_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	kom_sprintf_lite(str, "%s\n", rg_line);
	return 0;

err_set_rg:
	free(rg_line);
	return -1;
}

int mb_fmt_sam_hdr(kstring_t *str, const l2b_t *idx, const char *rg, const char *ver, int argc, char *argv[])
{
	int i, ret = 0;
	str->l = 0;
	kom_sprintf_lite(str, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
	if (idx)
		for (i = 0; i < idx->n_ctg; ++i)
			kom_sprintf_lite(str, "@SQ\tSN:%s\tLN:%ld\n", idx->ctg[i].name, idx->ctg[i].len);
	if (rg) ret = sam_write_rg_line(str, rg);
	kom_sprintf_lite(str, "@PG\tID:minibwa\tPN:minibwa");
	if (ver) kom_sprintf_lite(str, "\tVN:%s", ver);
	if (argc > 1) {
		kom_sprintf_lite(str, "\tCL:minibwa");
		for (i = 1; i < argc; ++i)
			kom_sprintf_lite(str, " %s", argv[i]);
	}
	return ret;
}

void mb_fmt_sam(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, int64_t opt_flag, int seg_idx)
{
}

void mb_format(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, int64_t opt_flag, int seg_idx)
{
	if (opt_flag & MB_F_SAM)
		mb_fmt_sam(km, s, l2b, t, n_seg, n_hit, hit, hit_idx, opt_flag, seg_idx);
	else
		mb_fmt_paf(s, l2b, t, &hit[seg_idx][hit_idx], opt_flag, n_seg, seg_idx);
}
