#include "mbpriv.h"
#include "kommon.h"

void mb_fmt_paf_basic(kstring_t *s, const l2b_t *l2b, int64_t qlen, const mb_hit_t *p, const char *qname)
{
	kom_sprintf_lite(s, "%s\t%ld\t%d\t%d\t%c\t%s\t%ld\t%ld\t%ld\t%d\t%d\t%d\ttp:A:%c\tcm:i:%d\n",
		qname, (long)qlen, p->qs, p->qe, p->rev? '-' : '+',
		l2b->ctg[p->tid].name, (long)l2b->ctg[p->tid].len, (long)p->ts, (long)p->te,
		p->mlen, p->blen, 255, p->parent == p->id? 'P' : 'S', p->cnt);
}
