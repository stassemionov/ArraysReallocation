#ifndef _REALLOCATION_H_
#define _REALLOCATION_H_

// This function makes inplace reallocating of standard-layout matrix to
// block representation by m x n blocks
double* block_reallocate_matrix(double* data_ptr,
    const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count);

// This function makes reallocating of standard-layout matrix to
// block representation by m x n blocks using a buffer array
double* get_reallocated(const double* data_ptr,
    const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count);

#endif  // _REALLOCATION_H_
