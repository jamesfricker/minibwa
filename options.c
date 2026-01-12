#include <string.h>
#include "minibwa.h"

void mb_mopt_init(mb_mopt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	// seeding options
	opt->min_len = 19;
	opt->max_sub_occ = 10;
	opt->max_occ = 500;
	// general algorithm options
	opt->bw = 500, opt->bw_long = 20000;
	opt->max_gap = 5000;
	// chaining options
	opt->max_chain_skip = 25;
	opt->max_chain_iter = 5000;
	opt->min_chain_score = 40;
	opt->rmq_inner_dist = 1000;
	opt->rmq_size_cap = 100000;
	opt->chn_pen_skip = 0.05f;
	// hit processing options
	opt->mask_level = 0.5f;
	opt->mask_len = 2147483647;
	opt->sub_diff = 0;
	opt->pri_ratio = 0.8f;
	opt->best_n = 5;
	// I/O options
	opt->n_thread = 4;
	opt->mb_size = 500000000;
}

int mb_set_preset(const char *preset, mb_mopt_t *opt)
{
	if (preset == 0) {
		mb_mopt_init(opt);
	} else if (strcmp(preset, "lr:hq") == 0) { // to be added
	} else if (strcmp(preset, "asm5") == 0) { // to be added
	} else if (strcmp(preset, "sr") == 0) {
		mb_mopt_init(opt);
		opt->flag |= MB_F_FRAG_MODE;
		opt->max_gap = 100;
		opt->bw = opt->bw_long = 100;
		opt->min_chain_score = 25;
		opt->max_occ = 5000;
		opt->mb_size = 50000000;
	} else {
		return -1;
	}
	return 0;
}
