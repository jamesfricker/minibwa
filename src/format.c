#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "mbpriv.h"
#include "kommon.h"

static char mb_rg_id[256];

static int str_contains_ci(const char *s, const char *needle);

static int str_eq_ci(const char *s, const char *t)
{
	int i;
	for (i = 0; s[i] && t[i]; ++i)
		if (tolower((unsigned char)s[i]) != tolower((unsigned char)t[i]))
			return 0;
	return s[i] == 0 && t[i] == 0;
}

static int str_starts_ci(const char *s, const char *prefix)
{
	int i;
	for (i = 0; prefix[i]; ++i)
		if (s[i] == 0 || tolower((unsigned char)s[i]) != tolower((unsigned char)prefix[i]))
			return 0;
	return 1;
}

static int is_primary_human_ctg(const char *name)
{
	const char *s = str_starts_ci(name, "chr")? name + 3 : name;
	char *end = 0;
	long n;
	if (str_eq_ci(s, "x") || str_eq_ci(s, "y") || str_eq_ci(s, "m") || str_eq_ci(s, "mt"))
		return 1;
	n = strtol(s, &end, 10);
	return end && *end == 0 && n >= 1 && n <= 22;
}

static const char *human_contig_class(const char *name)
{
	if (str_contains_ci(name, "hla")) return "hla";
	if (str_contains_ci(name, "decoy") || str_eq_ci(name, "hs37d5") || str_eq_ci(name, "hs38d1")) return "decoy";
	if (str_starts_ci(name, "chrUn") || str_contains_ci(name, "unplaced")) return "unplaced";
	if (str_contains_ci(name, "_random") || str_contains_ci(name, "random")) return "random";
	if (str_contains_ci(name, "_alt") || str_contains_ci(name, "alt_") || str_contains_ci(name, "alternate")) return "alt";
	if (is_primary_human_ctg(name)) return "primary";
	if (str_starts_ci(name, "GL") || str_starts_ci(name, "KI") || str_starts_ci(name, "JH") ||
		str_starts_ci(name, "KN") || str_starts_ci(name, "NT") || str_starts_ci(name, "NW"))
		return "unplaced";
	return "unknown";
}

static const char *human_region_class(const char *contig_class, double mappability)
{
	if (mappability < 0.5) return "low_mappability";
	if (strcmp(contig_class, "primary") == 0) return "high_confidence";
	if (strcmp(contig_class, "unknown") == 0) return "unknown";
	return "non_primary";
}

static void write_human_tags(kstring_t *s, const l2b_t *l2b, const mb_hit_t *r)
{
	const char *contig_class, *region_class;
	uint64_t aln_len, masked;
	double mappability = 1.0;
	int map_thousand = 1000;
	if (r == 0 || l2b == 0 || r->tid < 0 || r->tid >= (int64_t)l2b->n_ctg) return;
	contig_class = human_contig_class(l2b->ctg[r->tid].name);
	if (r->te > r->ts) {
		aln_len = r->te - r->ts;
		masked = l2b_mask_overlap(l2b, r->tid, r->ts, r->te);
		if (masked > aln_len) masked = aln_len;
		mappability = 1.0 - (double)masked / (double)aln_len;
		map_thousand = (int)(mappability * 1000.0 + 0.5);
		if (map_thousand < 0) map_thousand = 0;
		if (map_thousand > 1000) map_thousand = 1000;
	}
	region_class = human_region_class(contig_class, mappability);
	kom_sprintf_lite(s, "\tzc:Z:%s\tzm:f:%d.%d%d%d\tzh:Z:%s",
					  contig_class, map_thousand / 1000,
					  map_thousand / 100 % 10, map_thousand / 10 % 10, map_thousand % 10,
					  region_class);
}

/**************
 * PAF output *
 **************/

static inline void write_tags(kstring_t *s, const mb_hit_t *p)
{
	int32_t nm = p->blen - p->mlen + p->p->n_ambi;
	kom_sprintf_lite(s, "\tNM:i:%d\tAS:i:%d\tms:i:%d\tmd:i:%d", nm, p->p->dp_score, p->p->dp_max0, p->p->dp_max - p->p->dp_max2);
}

static inline void write_numt_tags(kstring_t *s, const mb_hit_t *p)
{
	if (p->numt_ambig)
		kom_sprintf_lite(s, "\tng:Z:chrM-nuclear");
}

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
	if (p->par) kom_sprintf_lite(s, "\tpa:Z:%s", mb_par_name(p->par));
	if (p->unmap)
		kom_sprintf_lite(s, "\tur:Z:unmap\tud:i:%d", p->unmap_max_depth);
	if (p->parent == p->id) kom_sprintf_lite(s, "\ts2:i:%d", p->subsc >= 0? p->subsc : 0);
	if (p->sv_blacklist) kom_sprintf_lite(s, "\tsb:Z:HMF_SV_BLACKLIST");
	if (p->problematic) kom_sprintf_lite(s, "\tgm:Z:GRC");
	if (opt_flag & MB_F_HUMAN_TAGS) write_human_tags(s, l2b, p);
	if (p->p) {
		write_tags(s, p);
		if (p->p->n_cigar > 0) {
			int32_t i;
			kom_sprintf_lite(s, "\tcg:Z:");
			for (i = 0; i < p->p->n_cigar; ++i)
				kom_sprintf_lite(s, "%d%c", p->p->cigar[i]>>4, MB_CIGAR_STR[p->p->cigar[i]&0xf]);
		}
		if (p->p->cs) kom_sprintf_lite(s, "\t%s", (char*)&p->p->cigar[p->p->n_cigar]);
	}
	write_numt_tags(s, p);
	if ((opt_flag & MB_F_COPY_COMMENT) && t->comment)
		kom_sprintf_lite(s, "\t%s", t->comment);
	kom_sprintf_lite(s, "\n");
}

/**************
 * SAM header *
 **************/

char *mb_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
			else if (*p == '\0') break;
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
	free(rg_line);
	return 0;

err_set_rg:
	free(rg_line);
	return -1;
}

int mb_fmt_sam_hdr(kstring_t *str, const l2b_t *idx, const char *rg, const char *ver, int argc, char *argv[], uint64_t opt_flag)
{
	int i, ret = 0;
	str->l = 0;
	kom_sprintf_lite(str, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
	if (idx)
		for (i = 0; i < idx->n_ctg; ++i)
			kom_sprintf_lite(str, "@SQ\tSN:%s\tLN:%ld\n", idx->ctg[i].name, idx->ctg[i].len);
	if (opt_flag & MB_F_HUMAN_TAGS) {
		kom_sprintf_lite(str, "@CO\tminibwa human-tags: zc=contig class from reference name; zm=unmasked alignment fraction from soft-masked reference; zh=human confidence class\n");
		kom_sprintf_lite(str, "@CO\tminibwa human-tags classes: zc primary,alt,hla,decoy,random,unplaced,unknown; zh high_confidence,low_mappability,non_primary,unknown\n");
	}
	if (rg) ret = sam_write_rg_line(str, rg);
	kom_sprintf_lite(str, "@PG\tID:minibwa\tPN:minibwa");
	if (ver) kom_sprintf_lite(str, "\tVN:%s", ver);
	if (argc > 1) {
		kom_sprintf_lite(str, "\tCL:minibwa");
		for (i = 0; i < argc; ++i)
			kom_sprintf_lite(str, " %s", argv[i]);
	}
	kom_sprintf_lite(str, "\n");
	return ret;
}

/**************
 * SAM output *
 **************/

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kom_roundup64(s->m);
		s->s = kom_realloc(char, s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

static void sam_write_sq(kstring_t *s, char *seq, int l, int rev, int comp)
{
	if (rev) {
		int i;
		str_enlarge(s, l);
		for (i = 0; i < l; ++i) {
			int c = seq[l - 1 - i];
			s->s[s->l + i] = c < 128 && comp? kom_comp_table[c] : c;
		}
		s->l += l;
	} else str_copy(s, seq, seq + l);
}

static inline const mb_hit_t *get_sam_pri(int n_hit, const mb_hit_t *hit)
{
	int i;
	for (i = 0; i < n_hit; ++i)
		if (hit[i].sam_pri)
			return &hit[i];
	assert(n_hit == 0);
	return NULL;
}

enum {
	MB_CTG_OTHER = 0,
	MB_CTG_PRIMARY = 1,
	MB_CTG_ALT = 2,
	MB_CTG_HLA = 3
};

static int ascii_lower(int c)
{
	return c >= 'A' && c <= 'Z'? c + ('a' - 'A') : c;
}

static int str_contains_ci(const char *s, const char *needle)
{
	int i;
	if (s == 0 || needle == 0 || needle[0] == 0) return 0;
	for (; *s; ++s) {
		for (i = 0; needle[i]; ++i)
			if (s[i] == 0 || ascii_lower((unsigned char)s[i]) != ascii_lower((unsigned char)needle[i]))
				break;
		if (needle[i] == 0) return 1;
	}
	return 0;
}

static int is_human_primary_name(const char *name)
{
	const char *p = name;
	int n = 0;
	if (p == 0) return 0;
	if (ascii_lower((unsigned char)p[0]) == 'c' &&
		ascii_lower((unsigned char)p[1]) == 'h' &&
		ascii_lower((unsigned char)p[2]) == 'r')
		p += 3;
	if (p[0] >= '1' && p[0] <= '9') {
		for (; *p >= '0' && *p <= '9'; ++p)
			n = n * 10 + (*p - '0');
		return *p == 0 && n >= 1 && n <= 22;
	}
	if ((p[0] == 'X' || p[0] == 'x' || p[0] == 'Y' || p[0] == 'y' || p[0] == 'M' || p[0] == 'm') && p[1] == 0)
		return 1;
	if ((p[0] == 'M' || p[0] == 'm') && (p[1] == 'T' || p[1] == 't') && p[2] == 0)
		return 1;
	return 0;
}

static int human_ctg_class(const l2b_t *l2b, int64_t tid)
{
	const char *name, *comm;
	if (l2b == 0 || tid < 0 || tid >= (int64_t)l2b->n_ctg) return MB_CTG_OTHER;
	name = l2b->ctg[tid].name;
	comm = l2b->ctg[tid].comm;
	if (str_contains_ci(name, "hla") || str_contains_ci(comm, "hla") ||
		str_contains_ci(name, "mhc") || str_contains_ci(comm, "mhc"))
		return MB_CTG_HLA;
	if (str_contains_ci(name, "_alt") || str_contains_ci(comm, "alternate locus") ||
		str_contains_ci(comm, "alt-scaffold") || str_contains_ci(comm, "alt scaffold"))
		return MB_CTG_ALT;
	return is_human_primary_name(name)? MB_CTG_PRIMARY : MB_CTG_OTHER;
}

static int same_query_locus(const mb_hit_t *a, const mb_hit_t *b)
{
	int st = a->qs > b->qs? a->qs : b->qs;
	int en = a->qe < b->qe? a->qe : b->qe;
	int al = a->qe - a->qs, bl = b->qe - b->qs;
	int min_l = al < bl? al : bl;
	return min_l > 0 && en > st && (en - st) * 5 >= min_l * 4;
}

static int human_alt_competitor(const l2b_t *l2b, const mb_hit_t *a, const mb_hit_t *b)
{
	int ac, bc;
	if (a == b || b->p == 0 || !same_query_locus(a, b)) return 0;
	ac = human_ctg_class(l2b, a->tid);
	bc = human_ctg_class(l2b, b->tid);
	if ((ac == MB_CTG_PRIMARY || ac == MB_CTG_ALT || ac == MB_CTG_HLA) &&
		(bc == MB_CTG_PRIMARY || bc == MB_CTG_ALT || bc == MB_CTG_HLA) &&
		(ac != MB_CTG_PRIMARY || bc != MB_CTG_PRIMARY))
		return 1;
	return 0;
}

static int sam_xa_candidate(const l2b_t *l2b, const mb_hit_t *r, const mb_hit_t *q, int r_i, int q_i, const mb_opt_t *opt)
{
	if (q_i == r_i || q->p == 0 || q->p->dp_max < (double)opt->xa_ratio * r->p->dp_max)
		return 0;
	if (q->parent == r_i) return 1;
	return (opt->flag & MB_F_HUMAN_ALT) && human_alt_competitor(l2b, r, q);
}

static void write_sam_cigar(kstring_t *s, int sam_flag, int in_tag, int qlen, const mb_hit_t *r, int64_t opt_flag)
{
	if (r->p == 0) {
		kom_sprintf_lite(s, "*");
	} else {
		uint32_t k, clip_len[2];
		clip_len[0] = r->rev? qlen - r->qe : r->qs;
		clip_len[1] = r->rev? r->qs : qlen - r->qe;
		if (in_tag) {
			int clip_char = (((sam_flag&0x800) || ((sam_flag&0x100) && (opt_flag&MB_F_2ND_SEQ))) &&
							 !(opt_flag&MB_F_SUPP_SOFT)) ? 5 : 4;
			kom_sprintf_lite(s, "\tCG:B:I");
			if (clip_len[0]) kom_sprintf_lite(s, ",%u", clip_len[0]<<4|clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				kom_sprintf_lite(s, ",%u", r->p->cigar[k]);
			if (clip_len[1]) kom_sprintf_lite(s, ",%u", clip_len[1]<<4|clip_char);
		} else {
			int clip_char = (((sam_flag&0x800) || ((sam_flag&0x100) && (opt_flag&MB_F_2ND_SEQ))) &&
							 !(opt_flag&MB_F_SUPP_SOFT)) ? 'H' : 'S';
			assert(clip_len[0] < qlen && clip_len[1] < qlen);
			if (clip_len[0]) kom_sprintf_lite(s, "%d%c", clip_len[0], clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				kom_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, MB_CIGAR_STR[r->p->cigar[k]&0xf]);
			if (clip_len[1]) kom_sprintf_lite(s, "%d%c", clip_len[1], clip_char);
		}
	}
}

void mb_fmt_sam(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, const mb_opt_t *opt, int seg_idx, int32_t mate_qlen)
{
	int flag, n_h = n_hit[seg_idx];
	int this_tid = -1, this_pos = -1;
	const mb_hit_t *h = hit[seg_idx], *r_prev = NULL, *r_next;
	const mb_hit_t *r = n_h > 0 && hit_idx < n_h && hit_idx >= 0? &h[hit_idx] : NULL;

	assert(n_seg == 1 || n_seg == 2);

	// find the primary of the previous and the next segments, if they are mapped
	if (n_seg > 1) {
		int next_sid = (seg_idx + 1) % n_seg;
		r_prev = r_next = get_sam_pri(n_hit[next_sid], hit[next_sid]);
	} else r_prev = r_next = NULL;

	// write QNAME and FLAG
	kom_sprintf_lite(s, "%s", t->name);
	flag = n_seg > 1? 0x1 : 0x0;
	if (r == 0) {
		flag |= 0x4;
	} else {
		if (r->rev) flag |= 0x10;
		if (r->parent != r->id) flag |= 0x100;
		else if (!r->sam_pri) flag |= 0x800;
	}
	if (n_seg > 1) {
		if (r && r->proper_pair) flag |= 0x2;
		if (seg_idx == 0) flag |= 0x40;
		else if (seg_idx == n_seg - 1) flag |= 0x80;
		if (r_next == NULL) flag |= 0x8;
		else if (r_next->rev) flag |= 0x20;
	}
	kom_sprintf_lite(s, "\t%d", flag);

	// write coordinate, MAPQ and CIGAR
	if (r == 0) {
		if (r_prev) {
			this_tid = r_prev->tid, this_pos = r_prev->ts;
			kom_sprintf_lite(s, "\t%s\t%d\t0\t*", l2b->ctg[this_tid].name, this_pos+1);
		} else kom_sprintf_lite(s, "\t*\t0\t0\t*");
	} else {
		this_tid = r->tid, this_pos = r->ts;
		kom_sprintf_lite(s, "\t%s\t%d\t%d\t", l2b->ctg[r->tid].name, r->ts+1, r->mapq);
		write_sam_cigar(s, flag, 0, t->l_seq, r, opt->flag);
	}

	// write mate positions
	if (n_seg > 1) {
		int tlen = 0;
		if (this_tid >= 0 && r_next) {
			if (this_tid == r_next->tid) {
				if (r) {
					int this_pos5 = r->rev? r->te - 1 : this_pos;
					int next_pos5 = r_next->rev? r_next->te - 1 : r_next->ts;
					tlen = next_pos5 - this_pos5;
				}
				kom_sprintf_lite(s, "\t=\t");
			} else kom_sprintf_lite(s, "\t%s\t", l2b->ctg[r_next->tid].name);
			kom_sprintf_lite(s, "%d\t", r_next->ts + 1);
		} else if (r_next) { // && this_tid < 0
			kom_sprintf_lite(s, "\t%s\t%d\t", l2b->ctg[r_next->tid].name, r_next->ts + 1);
		} else if (this_tid >= 0) { // && r_next == NULL
			kom_sprintf_lite(s, "\t=\t%d\t", this_pos + 1); // next segment will take r's coordinate
		} else kom_sprintf_lite(s, "\t*\t0\t"); // neither has coordinates
		if (tlen > 0) ++tlen;
		else if (tlen < 0) --tlen;
		kom_sprintf_lite(s, "%d\t", tlen);
	} else kom_sprintf_lite(s, "\t*\t0\t0\t");

	// write SEQ and QUAL
	if (r == 0) {
		sam_write_sq(s, t->seq, t->l_seq, 0, 0);
		kom_sprintf_lite(s, "\t");
		if (t->qual) sam_write_sq(s, t->qual, t->l_seq, 0, 0);
		else kom_sprintf_lite(s, "*");
	} else {
		if ((flag & 0x900) == 0 || (opt->flag & MB_F_SUPP_SOFT)) {
			sam_write_sq(s, t->seq, t->l_seq, r->rev, r->rev);
			kom_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual, t->l_seq, r->rev, 0);
			else kom_sprintf_lite(s, "*");
		} else if ((flag & 0x100) && !(opt->flag & MB_F_2ND_SEQ)){
			kom_sprintf_lite(s, "*\t*");
		} else {
			sam_write_sq(s, t->seq + r->qs, r->qe - r->qs, r->rev, r->rev);
			kom_sprintf_lite(s, "\t");
			if (t->qual) sam_write_sq(s, t->qual + r->qs, r->qe - r->qs, r->rev, 0);
			else kom_sprintf_lite(s, "*");
		}
	}

	// write tags
	if (mb_rg_id[0]) kom_sprintf_lite(s, "\tRG:Z:%s", mb_rg_id);
	if (n_seg > 2) kom_sprintf_lite(s, "\tFI:i:%d", seg_idx);
	if (r) {
		write_tags(s, r);
		if (r->par) kom_sprintf_lite(s, "\tpa:Z:%s", mb_par_name(r->par));
		if (r->sv_blacklist) kom_sprintf_lite(s, "\tsb:Z:HMF_SV_BLACKLIST");
		if (r->problematic) kom_sprintf_lite(s, "\tgm:Z:GRC");
		if (opt->flag & MB_F_HUMAN_TAGS) write_human_tags(s, l2b, r);
		if (r->unmap)
			kom_sprintf_lite(s, "\tur:Z:unmap\tud:i:%d", r->unmap_max_depth);
		write_numt_tags(s, r);
		// MC:Z mate CIGAR and MQ:i mate MAPQ; r_next is the mate's primary (see above).
		if (n_seg > 1 && r_next && r_next->p && r_next->p->n_cigar > 0 && mate_qlen > 0) {
			kom_sprintf_lite(s, "\tMC:Z:");
			write_sam_cigar(s, 0, 0, mate_qlen, r_next, opt->flag);
			kom_sprintf_lite(s, "\tMQ:i:%d", r_next->mapq);
		}
		if (r->p->cs) kom_sprintf_lite(s, "\t%s", (char*)&r->p->cigar[r->p->n_cigar]);
		if (r->parent == r->id && r->p && n_h > 1 && h && r >= h && r - h < n_h) { // supplementary aln may exist
			int i, r_i = r - h, n_sa = 0; // n_sa: number of SA fields
			for (i = 0; i < n_h; ++i)
				if (i != r_i && h[i].parent == h[i].id && h[i].p &&
					!(opt->xa_max > 0 && (opt->flag & MB_F_HUMAN_ALT) && human_alt_competitor(l2b, r, &h[i])))
					++n_sa;
			if (n_sa > 0) {
				kom_sprintf_lite(s, "\tSA:Z:");
				for (i = 0; i < n_h; ++i) {
					const mb_hit_t *q = &h[i];
					int l_M, l_I = 0, l_D = 0, clip5 = 0, clip3 = 0;
					if (r == q || q->parent != q->id || q->p == 0) continue;
					if (opt->xa_max > 0 && (opt->flag & MB_F_HUMAN_ALT) && human_alt_competitor(l2b, r, q)) continue;
					if (q->qe - q->qs < q->te - q->ts) l_M = q->qe - q->qs, l_D = (q->te - q->ts) - l_M;
					else l_M = q->te - q->ts, l_I = (q->qe - q->qs) - l_M;
					clip5 = q->rev? t->l_seq - q->qe : q->qs;
					clip3 = q->rev? q->qs : t->l_seq - q->qe;
					kom_sprintf_lite(s, "%s,%d,%c,", l2b->ctg[q->tid].name, q->ts+1, "+-"[q->rev]);
					if (clip5) kom_sprintf_lite(s, "%dS", clip5);
					if (l_M) kom_sprintf_lite(s, "%dM", l_M);
					if (l_I) kom_sprintf_lite(s, "%dI", l_I);
					if (l_D) kom_sprintf_lite(s, "%dD", l_D);
					if (clip3) kom_sprintf_lite(s, "%dS", clip3);
					kom_sprintf_lite(s, ",%d,%d;", q->mapq, q->blen - q->mlen + q->p->n_ambi);
				}
			}
			if (opt->xa_max > 0) {
				int i, n_xa = 0;
				for (i = 0; i < n_h; ++i)
					if (sam_xa_candidate(l2b, r, &h[i], r_i, i, opt))
						++n_xa;
				if (n_xa > 0) kom_sprintf_lite(s, "\tn2:i:%d", n_xa);
				if (n_xa > 0 && n_xa <= opt->xa_max) {
					kom_sprintf_lite(s, "\tXA:Z:");
					for (i = 0; i < n_h; ++i) {
						const mb_hit_t *q = &h[i];
						if (sam_xa_candidate(l2b, r, q, r_i, i, opt)) {
							kom_sprintf_lite(s, "%s,%c%d,", l2b->ctg[q->tid].name, "+-"[q->rev], q->ts+1);
							write_sam_cigar(s, 0, 0, t->l_seq, q, opt->flag);
							kom_sprintf_lite(s, ",%d;", q->blen - q->mlen + q->p->n_ambi);
						}
					}
				}
			}
		}
	}

	if ((opt->flag & MB_F_COPY_COMMENT) && t->comment)
		kom_sprintf_lite(s, "\t%s", t->comment);
	kom_sprintf_lite(s, "\n");
	s->s[s->l] = 0; // we always have room for an extra byte
}

void mb_format(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, const mb_opt_t *opt, int seg_idx, int32_t mate_qlen)
{
	if (!(opt->flag & MB_F_PAF))
		mb_fmt_sam(km, s, l2b, t, n_seg, n_hit, hit, hit_idx, opt, seg_idx, mate_qlen);
	else
		mb_fmt_paf(s, l2b, t, hit_idx >= 0? &hit[seg_idx][hit_idx] : 0, opt->flag, n_seg, seg_idx);
}
