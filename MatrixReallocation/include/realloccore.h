#ifndef _REALLOCCORE_H_
#define _REALLOCCORE_H_

#include "taskdata.h"

// Inplace reallocating of standard-layout matrix to
// block representation with b1 x b2 blocks.
const BlockReallocationInfo* standard_to_block_layout_reallocation(
    double* data_ptr,
    const TaskClass& data);

// Inplace reallocating of standard-layout matrix to
// block representation with b1 x b2 blocks,
// when cycles distribution data is specified by 'realloc_info' parameter.
void standard_to_block_layout_reallocation(
    double* data_ptr,
    const BlockReallocationInfo& realloc_info);

// Maps standard-layout source matrix to
// destination matrix with [b1 x b2] block layout.
double* map_with_block_layout(
    double* dst_ptr,
    const double* src_ptr,
    const TaskClass& task_info);

// Maps standard-layout source matrix to
// destination matrix with [b1 x b2, d1 x d2] double block layout.
double* map_with_double_block_layout(
    double* dst_ptr,
    const double* src_ptr,
    const TaskClass& task_info);

// Inplace reallocating of standard-layout matrix to
// double block representation with b1 x b2 main blocks
// and db1 x db2 small blocks.
const DoubleBlockReallocationInfo* standard_to_double_block_layout_reallocation(
    double* data_ptr,
    const TaskClass& data);
// Inplace reallocating of standard-layout matrix to
// double block representation with b1 x b2 main blocks
// and db1 x db2 small blocks,
// when cycles distribution data is specified by 'realloc_info' parameter.
void standard_to_double_block_layout_reallocation(
    double* data_ptr,
    const DoubleBlockReallocationInfo& realloc_info);

// Inplace reallocating of block-layout matrix to
// standard representation (row-major).
// Requires reallocation informaton (parameters + SDR vector).
double* block_to_standard_layout_reallocation(double* data_ptr,
    const BlockReallocationInfo& realloc_info);

// Inplace reallocating of double-block-layout matrix to
// standard representation (row-major).
// Requires reallocation informaton (parameters + SDR vector) for every
// big and small blocks sizes multiplicity case.
double* double_block_to_standard_layout_reallocation(double* data_ptr,
    const DoubleBlockReallocationInfo& realloc_info);

#endif  // _REALLOCCORE_H_
