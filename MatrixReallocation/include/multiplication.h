#ifndef _MULTIPLICATION_H_
#define _MULTIPLICATION_H_

#ifdef __cplusplus
extern "C" {
#endif

void matrix_multiplication_double_block(double* gen_matrix,
    const double* left_matrix,
    const double* right_matrix,
    const int N1, const int N2, const int N3,
    const int B1, const int B2, const int B3,
    const int D1, const int D2, const int D3);

void matrix_multiplication_block(double* gen_matrix,
    const double* left_matrix,
    const double* right_matrix,
    const int N1, const int N2, const int N3,
    const int B1, const int B2, const int B3);

void matrix_multiplication_double_tiled(double* gen_matrix,
    const double* left_matrix,
    const double* right_matrix,
    const int N1, const int N2, const int N3,
    const int B1, const int B2, const int B3,
    const int D1, const int D2, const int D3);

void matrix_multiplication_tiled(double* gen_matrix,
    const double* left_matrix,
    const double* right_matrix,
    const int N1, const int N2, const int N3,
    const int B1, const int B2, const int B3);

void matrix_multiplication_standard(double* gen_matrix,
    const double* src_left_matrix,
    const double* src_right_matrix,
    const int N1, const int N2, const int N3);

#ifdef __cplusplus
    }
#endif

#endif  // _MULTIPLICATION_H_
