#ifndef _TESTING_H_
#define _TESTING_H_

#include "taskdata.h"

struct Blocks
{
    int main_block;
    int small_block;
};

// Tool for matrix multiplication performance testing.
// Performs tiled and double tiled multiplication with standard allocation
// and tiled and double tiled multiplication with block and
// double block allocations correspondingly.
void matrix_multiplication_tests(const TaskClass& params_left,
    const TaskClass& params_right,
    const bool       console_info_output = false);

void floyd_test(const TaskClass& parameters,
    const bool console_info_output);

void qralg_test(const TaskClass& parameters,
    const bool console_info_output);

void correctness_test(const TaskClass& params_left,
    const TaskClass& params_right,
    const bool       console_info_output = false);

void reallocation_test(const TaskClass& parameters,
                       const bool console_info_output = false);


// Tool for block size selection for square matrix block multiplication.
// Require 40*N*N bytes of additional memory.
int select_optimal_block_size_multiplication(
    const int N,
    const int start_block_val,
    const int end_block_val,
    const int selection_step);

// Tool for double block size selection for square matrix block multiplication.
// Requires 40*N*N bytes of additional memory.
// Returns structs that is pair - main block size and small block size.
Blocks select_optimal_double_block_size_multiplication(
    const int N,
    const int start_block_val,
    const int end_block_val,
    const int selection_block_step,
    const int start_double_block_val,
    const int end_double_block_val,
    const int selection_double_block_step,
    const bool         console_info_output = false);

int select_optimal_block_size_floyd(
    const int N,
    const int start_block_val,
    const int end_block_val,
    const int selection_step);

#endif  // _TESTING_H_
