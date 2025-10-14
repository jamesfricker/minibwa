#ifndef MB_SYS_H
#define MB_SYS_H

#ifdef __cplusplus
extern "C" {
#endif

double mb_cputime(void);
double mb_realtime(void);
long mb_peakrss(void);

#ifdef __cplusplus
}
#endif

#endif
