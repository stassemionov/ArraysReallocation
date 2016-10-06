#ifndef _QRALG_H_
#define _QRALG_H_

#include "taskdata.h"

// Standard QR-algorithm with WY-decomposition.
double* QR_WY_standard(double* A, const int N, const int r);

// Standard QR-algorithm with WY-decomposition.
// Matrix multiplications are produced with tiling.
double* QR_WY_tiled(double* A, const TaskData& parameters);

// Standard QR-algorithm with WY-decomposition.
// Matrix multiplications are produced with double tiling.
double* QR_WY_double_tiled(double* A, const TaskData& parameters);

// Standard QR-algorithm with WY-decomposition.
// Matrix A must has block layout.
double* QR_WY_block(double* A, const TaskData& parameters);

// Standard QR-algorithm with WY-decomposition.
// Matrix A must has double block layout.
double* QR_WY_double_block(double* A, const TaskData& parameters);

#endif  // _QRALG_H_
