#ifndef _FLOYDALG_H_
#define _FLOYDALG_H_

double* floyd_standard(double* A, const int N);

double* floyd_tiled(double* A, const int N, const int B);

double* floyd_block(double* data, const int N, const int B);

#endif  // _FLOYDALG_H_
