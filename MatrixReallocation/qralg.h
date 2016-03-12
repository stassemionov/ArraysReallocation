#ifndef _QRALG_H_
#define _QRALG_H_

#include "taskdata.h"

// Стандартный QR-алгоритм с WY-разложением.
double* QR_WY_standard(double* A, const int N, const int r);

// QR-алгоритм с WY-разложением.
// Умножения матриц производятся с тайлингом.
double* QR_WY_tiled(double* A, const TaskData& parameters);

// QR-алгоритм с WY-разложением.
// Матрица А предполагается размещенной блочно.
double* QR_WY_block(double* A, const TaskData& parameters);

#endif  // _QRALG_H_