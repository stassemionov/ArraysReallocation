#ifndef _REALLOCATION_H_
#define _REALLOCATION_H_

#include "taskdata.h"

// Inplace reallocating of standard-layout matrix to
// block representation with b1 x b2 blocks.
double* standard_to_block_layout_reallocation(double* data_ptr,
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
double* standard_to_double_block_layout_reallocation(double* data_ptr,
                                                     const TaskClass& data);

#endif  // _REALLOCATION_H_
