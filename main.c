#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwt.h"
#include "sys.h"
#include "utils.h"
#include "ketopt.h"

#define MB_VERSION "0.0"

int main_raw2bwt(int argc, char *argv[]);
int main_test(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minibwt <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  raw2bwt    recode bwtgen raw BWT\n");
	fprintf(fp, "  version    print the version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	mb_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "raw2bwt") == 0) ret = main_raw2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "test") == 0) ret = main_test(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("%s\n", MB_VERSION);
		return 0;
	} else {
		fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
		return 1;
	}

	if (mb_verbose >= 3 && argc > 2 && ret == 0) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MB_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, mb_realtime(), mb_cputime(), mb_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
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

int main_test(int argc, char *argv[])
{
	uint64_t i, n = 10000, x = 11;
	mb_bwt_t *bwt;
	if (argc == 1) {
		fprintf(stderr, "Usage: minibwa test <in.bwt>\n");
		return 1;
	}
	if (argc >= 3) n = atol(argv[2]);

	bwt = mb_bwt_load_raw(argv[1]);
	if (n > 0) {
		double t = mb_cputime();
		for (i = 0; i < n; ++i) {
			uint64_t k = mb_splitmix64(&x) % bwt->seq_len;
			#if 1
			uint64_t cnt[4];
			mb_bwt_rank1a(bwt, k, cnt);
			printf("%lld\n", cnt[1]);
			#else
			printf("%lld\n", mb_bwt_rank11(bwt, k, 1));
			#endif
		}
		fprintf(stderr, "t = %.3f\n", mb_cputime() - t);
	} else {
		uint64_t k = 10, cnt[4];
		int c;
		mb_bwt_rank1a(bwt, k, cnt);
		for (c = 0; c < 4; ++c)
			printf("%lld\n", cnt[c]);
	}
	mb_bwt_destroy(bwt);
	return 0;
}
