﻿#ifndef _REALLOCATION_H_
#define _REALLOCATION_H_

#include "taskdata.h"

// Inplace reallocating of standard-layout matrix to
// block representation with b1 x b2 blocks.
BlockReallocationInfo* standard_to_block_layout_reallocation(
                                                double* data_ptr,
                                                const TaskClass& data);

// Reallocating of standard-layout matrix to
// block representation with b1 x b2 blocks using a buffer array.
double* standard_to_block_layout_reallocation_buf(const double* data_ptr,
                                                  const TaskClass& data);

// Reallocating of standard-layout matrix to
// double block representation with b1 x b2 main blocks
// and d1 x d2 subblocks, using a buffer array.
double* standard_to_double_block_layout_reallocation_buf(
    const double* data_ptr, const TaskClass& data);

// Inplace reallocating of standard-layout matrix to
// double block representation with b1 x b2 main blocks
// and db1 x db2 small blocks.
DoubleBlockReallocationInfo* standard_to_double_block_layout_reallocation(
                                                     double* data_ptr,
                                                     const TaskClass& data);

// Inplace reallocating of block-layout matrix to
// standard representation (row-major).
// Requires reallocation informaton (parameters + SDR vector).
double* block_to_standard_layout_reallocation(double* data_ptr,
    const BlockReallocationInfo* realloc_info);

// Inplace reallocating of double-block-layout matrix to
// standard representation (row-major).
// Requires reallocation informaton (parameters + SDR vector) for every
// big and small blocks sizes multiplicity case.
double* double_block_to_standard_layout_reallocation(double* data_ptr,
    const DoubleBlockReallocationInfo* realloc_info);

#endif  // _REALLOCATION_H_
