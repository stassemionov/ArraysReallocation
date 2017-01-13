#ifndef _REALLOCATION_H_
#define _REALLOCATION_H_

// Must be called before first using of reallocation functions
void InitDispatchSystem();

// Must be called after last using of reallocation functions
void TurnOffDispatchSystem();

double* standard_to_block_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2);

double* standard_to_double_block_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2);

double* block_to_standard_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2);

double* double_block_to_standard_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2);

double* standard_to_transposed_double_block_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2);

double* transposed_double_block_to_standard_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2);

#endif  // _REALLOCATION_H_
