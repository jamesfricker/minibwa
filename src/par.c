#include <string.h>
#include "mbpriv.h"

typedef struct {
	int32_t id;
	int64_t x_st, x_en;
	int64_t y_st, y_en;
} mb_par_intv_t;

static const mb_par_intv_t mb_par_grch37[] = {
	{ 1, 60000, 2699520, 10000, 2649520 },
	{ 2, 154931043, 155260560, 59034049, 59363566 }
};

static const mb_par_intv_t mb_par_grch38[] = {
	{ 1, 10000, 2781479, 10000, 2781479 },
	{ 2, 155701382, 156030895, 56887902, 57217415 }
};

static const mb_par_intv_t *mb_par_table(int32_t asm_id, int32_t *n)
{
	if (asm_id == 37) {
		*n = sizeof(mb_par_grch37) / sizeof(mb_par_grch37[0]);
		return mb_par_grch37;
	} else if (asm_id == 38) {
		*n = sizeof(mb_par_grch38) / sizeof(mb_par_grch38[0]);
		return mb_par_grch38;
	}
	*n = 0;
	return 0;
}

static int mb_is_x_name(const char *s)
{
	return strcmp(s, "chrX") == 0 || strcmp(s, "X") == 0 ||
		   strcmp(s, "NC_000023.10") == 0 || strcmp(s, "NC_000023.11") == 0;
}

static int mb_is_y_name(const char *s)
{
	return strcmp(s, "chrY") == 0 || strcmp(s, "Y") == 0 ||
		   strcmp(s, "NC_000024.9") == 0 || strcmp(s, "NC_000024.10") == 0;
}

void mb_par_init(l2b_t *l2b)
{
	int64_t i, x = -1, y = -1;
	if (l2b == 0) return;
	l2b->par_asm = -1;
	l2b->par_x_tid = l2b->par_y_tid = -1;
	for (i = 0; i < l2b->n_ctg; ++i) {
		if (mb_is_x_name(l2b->ctg[i].name)) x = i;
		else if (mb_is_y_name(l2b->ctg[i].name)) y = i;
	}
	if (x < 0 || y < 0) return;
	if (l2b->ctg[x].len == 155270560 && l2b->ctg[y].len == 59373566)
		l2b->par_asm = 37;
	else if (l2b->ctg[x].len == 156040895 && l2b->ctg[y].len == 57227415)
		l2b->par_asm = 38;
	else return;
	l2b->par_x_tid = (int32_t)x;
	l2b->par_y_tid = (int32_t)y;
}

static int mb_par_project(const l2b_t *l2b, const mb_hit_t *h, int64_t *ps, int64_t *pe)
{
	const mb_par_intv_t *t;
	int32_t i, n;
	if (l2b == 0 || h == 0 || l2b->par_asm < 0) return 0;
	t = mb_par_table(l2b->par_asm, &n);
	for (i = 0; i < n; ++i) {
		int64_t st = -1, en = -1;
		if (h->tid == l2b->par_x_tid) st = t[i].x_st, en = t[i].x_en;
		else if (h->tid == l2b->par_y_tid) st = t[i].y_st, en = t[i].y_en;
		if (st >= 0 && h->ts >= st && h->te <= en) {
			*ps = h->ts - st;
			*pe = h->te - st;
			return t[i].id;
		}
	}
	return 0;
}

void mb_mark_par_hits(const l2b_t *l2b, int32_t n, mb_hit_t *h)
{
	int32_t i;
	for (i = 0; i < n; ++i) {
		int64_t ps, pe;
		h[i].par = mb_par_project(l2b, &h[i], &ps, &pe);
	}
}

int mb_par_equiv(const l2b_t *l2b, const mb_hit_t *a, const mb_hit_t *b)
{
	int64_t as, ae, bs, be;
	if (a == 0 || b == 0 || a->par == 0 || a->par != b->par) return 0;
	if (a->tid == b->tid || a->rev != b->rev || a->qs != b->qs || a->qe != b->qe) return 0;
	if (mb_par_project(l2b, a, &as, &ae) != a->par) return 0;
	if (mb_par_project(l2b, b, &bs, &be) != b->par) return 0;
	return as == bs && ae == be;
}

static void mb_par_recalc_parent_evidence(const l2b_t *l2b, int32_t n, mb_hit_t *h, int32_t sub_diff)
{
	int32_t i;
	for (i = 0; i < n; ++i) {
		if (h[i].parent == h[i].id) {
			h[i].subsc = 0;
			h[i].n_sub = 0;
		}
	}
	for (i = 0; i < n; ++i) {
		mb_hit_t *c = &h[i], *p;
		if (c->parent < 0 || c->parent == c->id || c->parent >= n) continue;
		p = &h[c->parent];
		if (p->parent != p->id || mb_par_equiv(l2b, c, p)) continue;
		if (p->subsc < c->score) p->subsc = c->score;
		if (p->p && c->p) {
			if (p->p->dp_max2 < c->p->dp_max) p->p->dp_max2 = c->p->dp_max;
			if (p->p->dp_max - c->p->dp_max <= sub_diff) ++p->n_sub;
		}
	}
}

void mb_par_resolve(const l2b_t *l2b, int32_t n, mb_hit_t *h, int32_t sub_diff)
{
	int32_t i, j, changed = 0;
	if (l2b == 0 || l2b->par_asm < 0) return;
	for (i = 0; i < n; ++i) {
		for (j = i + 1; j < n; ++j) {
			int32_t k, x, y;
			if (!mb_par_equiv(l2b, &h[i], &h[j])) continue;
			if (h[i].tid == l2b->par_x_tid && h[j].tid == l2b->par_y_tid) x = i, y = j;
			else if (h[j].tid == l2b->par_x_tid && h[i].tid == l2b->par_y_tid) x = j, y = i;
			else continue;
			if (h[x].parent != h[x].id && h[y].parent != h[y].id && h[x].parent != h[y].id && h[y].parent != h[x].id)
				continue;
			for (k = 0; k < n; ++k)
				if (h[k].parent == h[y].id)
					h[k].parent = h[x].id;
			h[x].parent = h[x].id;
			h[y].parent = h[x].id;
			changed = 1;
		}
	}
	if (changed) mb_par_recalc_parent_evidence(l2b, n, h, sub_diff);
}

const char *mb_par_name(int32_t par)
{
	return par == 1? "PAR1" : par == 2? "PAR2" : 0;
}
