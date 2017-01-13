#ifndef _MAPPINGFUNCS_H_
#define _MAPPINGFUNCS_H_

#include "taskdata.h"

double* map_with_block_layout(
    double* dst_ptr,
    const double* src_ptr,
    const int N1, const int N2,
    const int B1, const int B2);

double* map_with_double_block_layout(
    double* dst_ptr,
    const double* src_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2);

double* map_with_transposed_double_block_layout(
    double* dst_ptr,
    const double* src_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2);

#endif  // _MAPPINGFUNCS_H_
