#ifndef _FLOYDALG_H_
#define _FLOYDALG_H_

#ifdef __cplusplus
extern "C" {
#endif

double* floyd_standard(double* A, const int N);

double* floyd_tiled(double* A, const int N, const int B);

double* floyd_block(double* data, const int N, const int B);

#ifdef __cplusplus
    }
#endif

#endif  // _FLOYDALG_H_
