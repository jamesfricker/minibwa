#include <stdlib.h>
#include <stdio.h>
#include "bwt.h"
#include "l2bit.h"
#include "libsais64.h"
#include "kommon.h"
#include "ketopt.h"

mb_bwt_t *mb_bwt_libsais(const l2b_t *l2b, int both_strand, int n_thread)
{
	const int fs = 10000;
	uint8_t *seq;
	int64_t i, j, *a, primary, len;
	mb_bwt_t *bwt;

	len = both_strand? l2b->tot_len * 2 : l2b->tot_len;
	seq = kom_malloc(uint8_t, len);
	a = kom_malloc(int64_t, len + fs + 1);
	for (i = 0, j = 0; i < l2b->tot_len; ++i, ++j)
		seq[j] = l2b_get0(l2b, i);
	if (both_strand)
		for (i = l2b->tot_len - 1; i >= 0; --i, ++j)
			seq[j] = 3 - l2b_get0(l2b, i);
#ifdef LIBSAIS_OPENMP
    primary = libsais64_bwt_omp(seq, seq, a, len, fs, 0, n_thread);
#else
    primary = libsais64_bwt(seq, seq, a, len, fs, 0);
#endif
	free(a);
	bwt = mb_bwt_init_from_raw(1, seq, len, primary);
	free(seq);
	return bwt;
}

int main_fa2bit(int argc, char *argv[])
{
	l2b_t *l2b;
	int out_pac = 0, both_strand = 0;
	uint64_t seed = 11;
	ketopt_t o = KETOPT_INIT;
	int c;
	while ((c = ketopt(&o, argc, argv, 1, "s:p2", 0)) >= 0) {
		if (c == 's') seed = atol(o.arg);
		else if (c == 'p') out_pac = 1;
		else if (c == '2') both_strand = 1;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minibwa fa2bit [options] <in.fa> <out.l2b>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -s INT    random seed [%lu]\n", (unsigned long)seed);
		fprintf(stderr, "  -p        output the BWA pac format\n");
		fprintf(stderr, "  -2        output both strands (effective with -p)\n");
		return 1;
	}
	l2b = l2b_import(argv[o.ind], seed);
	if (out_pac)
		l2b_save_pac(argv[o.ind+1], l2b, both_strand);
	else
		l2b_save(argv[o.ind+1], l2b);
	l2b_destroy(l2b);
	return 0;
}

int main_genraw(int argc, char *argv[])
{
	extern void mb_bwtgen(const char *fn_pac, const char *fn_bwt, int block_size);
	ketopt_t o = KETOPT_INIT;
	int c, block_size = 10000000;
	while ((c = ketopt(&o, argc, argv, 1, "b:", 0)) >= 0) {
		if (c == 'b') block_size = kom_parse_num(o.arg, 0);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minibwa genraw [options] <in.pac> <out.raw-bwt>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b NUM      block size [10m]\n");
		return 1;
	}
	mb_bwtgen(argv[o.ind], argv[o.ind+1], block_size);
	return 0;
}

int main_raw2bwt(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	if (argc < 3) {
		printf("Usage: minibwa raw2bwt <raw.bwt> <recode.bwt>\n");
		return 1;
	}
	bwt = mb_bwt_load_raw(argv[1]);
	mb_bwt_save(argv[2], bwt);
	mb_bwt_destroy(bwt);
	return 0;
}

int main_genbwt(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, n_thread = 1, both_strand = 1;
	mb_bwt_t *bwt;
	l2b_t *l2b;
	while ((c = ketopt(&o, argc, argv, 1, "1t:", 0)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		else if (c == '1') both_strand = 0;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minibwa genbwt [options] <in.l2b> <out.bwt>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -1          forward strand only\n");
#ifdef LIBSAIS_OPENMP
		fprintf(stderr, "  -t INT      number of threads\n");
#endif
		return 1;
	}
	l2b = l2b_load(argv[o.ind]);
	kom_assert(l2b, "failed to open the input file.");
	bwt = mb_bwt_libsais(l2b, both_strand, n_thread);
	l2b_destroy(l2b);
	mb_bwt_save(argv[o.ind+1], bwt);
	mb_bwt_destroy(bwt);
	return 0;
}

int main_gensa(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	int c, sa_bit = 5, is_raw = 0;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "rb:", 0)) >= 0) {
		if (c == 'b') sa_bit = atoi(o.arg);
		else if (c == 'r') is_raw = 1;
	}
	if (argc < 3) {
		fprintf(stderr, "Usage: minibwa gensa [options] <in.bwt> <out.bwt>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b INT    sample rate at 1/(1<<INT) [%d]\n", sa_bit);
		fprintf(stderr, "  -r        input BWT in the raw BWA format\n");
		return 1;
	}
	bwt = is_raw? mb_bwt_load_raw(argv[o.ind]) : mb_bwt_load(argv[o.ind]);
	mb_bwt_gen_sa(bwt, sa_bit);
	mb_bwt_save(argv[o.ind+1], bwt);
	mb_bwt_destroy(bwt);
	return 0;
}
