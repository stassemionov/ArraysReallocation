#ifndef _TESTING_H_
#define _TESTING_H_

struct Blocks
{
    unsigned int main_block;
    unsigned int small_block;
};

// Tool for block size selection for square matrix block multiplication.
// Require 40*N*N bytes of additional memory.
unsigned int select_optimal_block_size(const unsigned int N,
    const unsigned int start_block_val,
    const unsigned int end_block_val,
    const unsigned int selection_step);

// Tool for double block size selection for square matrix block multiplication.
// Require 40*N*N bytes of additional memory.
// Returns structs that is pair - main block size and small block size.
Blocks select_optimal_double_block_size(const unsigned int N,
    const unsigned int start_block_val,
    const unsigned int end_block_val,
    const unsigned int selection_block_step,
    const unsigned int start_double_block_val,
    const unsigned int end_double_block_val,
    const unsigned int selection_double_block_step,
    const bool         console_info_output = false);

#endif  // _TESTING_H_
