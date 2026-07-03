#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "mbpriv.h"
#include "ksort.h"

#define key_unmap(a) ((uint64_t)(a).st)
KRADIX_SORT_INIT(unmap, mb_unmap_intv_t, key_unmap, 8)

static int split_fields(char *line, char **field, int max_field)
{
	int n = 0;
	char *p = line;
	while (*p) {
		while (*p == '\t' || *p == ',' || *p == ' ' || *p == '\r' || *p == '\n') ++p;
		if (*p == 0) break;
		if (n == max_field) break;
		field[n++] = p;
		while (*p && *p != '\t' && *p != ',' && *p != ' ' && *p != '\r' && *p != '\n') ++p;
		if (*p == 0) break;
		*p++ = 0;
	}
	return n;
}

static int field_index(char **field, int n, const char *name)
{
	int i;
	for (i = 0; i < n; ++i)
		if (strcmp(field[i], name) == 0)
			return i;
	return -1;
}

static int ctg_name_match(const char *a, const char *b)
{
	if (strcmp(a, b) == 0) return 1;
	if (strncmp(a, "chr", 3) == 0 && strcmp(a + 3, b) == 0) return 1;
	if (strncmp(b, "chr", 3) == 0 && strcmp(a, b + 3) == 0) return 1;
	return 0;
}

static int32_t find_ctg(const l2b_t *l2b, const char *name)
{
	int32_t i;
	for (i = 0; i < (int32_t)l2b->n_ctg; ++i)
		if (ctg_name_match(l2b->ctg[i].name, name))
			return i;
	return -1;
}

static void add_region(mb_unmap_regions_t *r, int32_t tid, int64_t st, int64_t en, int32_t max_depth)
{
	mb_unmap_ctg_t *ctg = &r->ctg[tid];
	mb_unmap_intv_t *iv;
	if (st < 0) st = 0;
	if (en <= st) return;
	kom_grow(mb_unmap_intv_t, ctg->a, ctg->n, ctg->m);
	iv = &ctg->a[ctg->n++];
	iv->st = st, iv->en = en, iv->max_en = en, iv->max_depth = max_depth;
	++r->n_regions;
}

mb_unmap_regions_t *mb_unmap_regions_load(const char *fn, const l2b_t *l2b)
{
	gzFile fp;
	char line[65536], *field[64];
	int n, chr_i, st_i, en_i, depth_i;
	mb_unmap_regions_t *r;

	fp = gzopen(fn, "rb");
	if (fp == 0) return 0;
	if (gzgets(fp, line, sizeof(line)) == 0) {
		gzclose(fp);
		return 0;
	}
	n = split_fields(line, field, 64);
	chr_i = field_index(field, n, "Chromosome");
	st_i = field_index(field, n, "PosStart");
	en_i = field_index(field, n, "PosEnd");
	depth_i = field_index(field, n, "MaxDepth");
	if (chr_i < 0 || st_i < 0 || en_i < 0 || depth_i < 0) {
		gzclose(fp);
		return 0;
	}

	r = kom_calloc(mb_unmap_regions_t, 1);
	r->n_ctg = (int32_t)l2b->n_ctg;
	r->ctg = kom_calloc(mb_unmap_ctg_t, r->n_ctg);
	r->source = kom_strdup(fn);
	while (gzgets(fp, line, sizeof(line)) != 0) {
		int32_t tid, max_depth;
		int64_t pos_start, pos_end;
		if (line[0] == 0 || line[0] == '#' || line[0] == '\n') continue;
		n = split_fields(line, field, 64);
		if (n <= depth_i) continue;
		tid = find_ctg(l2b, field[chr_i]);
		if (tid < 0) continue;
		pos_start = strtoll(field[st_i], 0, 10);
		pos_end = strtoll(field[en_i], 0, 10);
		max_depth = (int32_t)strtol(field[depth_i], 0, 10);
		// HMF PosStart/PosEnd are 1-based inclusive; minibwa hits are 0-based half-open.
		add_region(r, tid, pos_start - 1, pos_end, max_depth);
	}
	gzclose(fp);

	for (n = 0; n < r->n_ctg; ++n) {
		mb_unmap_ctg_t *ctg = &r->ctg[n];
		int32_t i;
		if (ctg->n > 1)
			radix_sort_unmap(ctg->a, ctg->a + ctg->n);
		for (i = 0; i < ctg->n; ++i)
			ctg->a[i].max_en = i == 0 ? ctg->a[i].en
				: (ctg->a[i].en > ctg->a[i-1].max_en ? ctg->a[i].en : ctg->a[i-1].max_en);
	}
	return r;
}

void mb_unmap_regions_destroy(mb_unmap_regions_t *r)
{
	int32_t i;
	if (r == 0) return;
	for (i = 0; i < r->n_ctg; ++i)
		free(r->ctg[i].a);
	free(r->ctg);
	free(r->source);
	free(r);
}

int32_t mb_unmap_hit_max_depth(const mb_unmap_regions_t *r, const mb_hit_t *h)
{
	const mb_unmap_ctg_t *ctg;
	int32_t lo, hi, mid, end, i, max_depth = 0;
	if (r == 0 || h == 0 || h->tid < 0 || h->tid >= r->n_ctg || h->te <= h->ts) return 0;
	ctg = &r->ctg[h->tid];
	if (ctg->n == 0) return 0;
	lo = 0, hi = ctg->n;
	while (lo < hi) {
		mid = (lo + hi) >> 1;
		if (ctg->a[mid].st < h->te) lo = mid + 1;
		else hi = mid;
	}
	end = lo;
	lo = 0, hi = end;
	while (lo < hi) {
		mid = (lo + hi) >> 1;
		if (ctg->a[mid].max_en > h->ts) hi = mid;
		else lo = mid + 1;
	}
	for (i = lo; i < end; ++i)
		if (ctg->a[i].en > h->ts && max_depth < ctg->a[i].max_depth)
			max_depth = ctg->a[i].max_depth;
	return max_depth;
}

void mb_apply_unmap_regions(const mb_unmap_regions_t *r, int32_t n_hit, mb_hit_t *hit)
{
	int32_t i;
	if (r == 0 || n_hit <= 0 || hit == 0) return;
	for (i = 0; i < n_hit; ++i) {
		hit[i].unmap_max_depth = mb_unmap_hit_max_depth(r, &hit[i]);
		hit[i].unmap = hit[i].unmap_max_depth > 0;
	}
}
