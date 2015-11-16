#ifndef _REALLOCATING_H_
#define _REALLOCATING_H_

// This function makes inplace reallocating of standard-layout matrix to
// block representation by m x n blocks
double* block_reallocate_matrix(double* A, int M, int N, int m, int n);

// This function makes reallocating of standard-layout matrix to
// block representation by m x n blocks using a buffer array
double* buffer_reallocating(double* A, int M, int N, int m, int n);

#endif _REALLOCATING_H_