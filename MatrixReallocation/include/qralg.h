#ifndef _QRALG_H_
#define _QRALG_H_

#ifdef __cplusplus
extern "C" {
#endif

// Standard QR-algorithm with WY-decomposition.
double* QR_WY_standard(double* A, const int N, const int r);

// Standard QR-algorithm with WY-decomposition.
// Matrix multiplications are produced with tiling.
double* QR_WY_tiled(double* A, const int N, const int b1, const int b2);

// Standard QR-algorithm with WY-decomposition.
// Matrix multiplications are produced with double tiling.
double* QR_WY_double_tiled(double* A, 
    const int N, const int b1, const int b2, const int db1, const int db2);

// Standard QR-algorithm with WY-decomposition.
// Matrix A must has block layout.
double* QR_WY_block(double* A, const int N, const int b1, const int b2);

// Standard QR-algorithm with WY-decomposition.
// Matrix A must has double block layout.
double* QR_WY_double_block(double* A,
    const int N, const int b1, const int b2, const int db1, const int db2);

#ifdef __cplusplus
    }
#endif

#endif  // _QRALG_H_
