#ifndef HUMAN_RESOURCES_H
#define HUMAN_RESOURCES_H

#include "l2bit.h"

#ifdef __cplusplus
extern "C" {
#endif

int mb_human_resources_import(const char *src, const char *prefix, const l2b_t *l2b);
void mb_human_resources_warn_mismatch(const char *prefix, const l2b_t *l2b);

#ifdef __cplusplus
}
#endif

#endif
