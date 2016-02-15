#ifndef _MULTIPLICATION_H_
#define _MULTIPLICATION_H_

#include "taskdata.h"

void block_matrix_multiplication_double_tiled(double* gen_matrix,
    const double* left_matrix,
    const double* right_matrix,
    const TaskClass& left_mat_data,
    const TaskClass& right_mat_data);

void block_matrix_multiplication_tiled(double* gen_matrix,
    const double* left_matrix,
    const double* right_matrix,
    const TaskClass& left_mat_data,
    const TaskClass& right_mat_data);

void matrix_multiplication_double_tiled(double* gen_matrix,
    const double* left_matrix,
    const double* right_matrix,
    const TaskClass& left_mat_data,
    const TaskClass& right_mat_data);

void matrix_multiplication_tiled(double* gen_matrix,
    const double* left_matrix,
    const double* right_matrix,
    const TaskClass& left_mat_data,
    const TaskClass& right_mat_data);

void matrix_multiplication_standard(double* gen_matrix,
    const double* src_left_matrix,
    const double* src_right_matrix,
    const int N1, const int N2, const int N3);

#endif  // _MULTIPLICATION_H_
