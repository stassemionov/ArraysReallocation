#ifndef _REALLOCATION_H_
#define _REALLOCATION_H_

// Inplace reallocating of standard-layout matrix to
// block representation with b1 x b2 blocks.
double* standard_to_block_layout_reallocation(double* data_ptr,
    const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count);

// Reallocating of standard-layout matrix to
// block representation with b1 x b2 blocks using a buffer array.
double* standard_to_block_layout_reallocation_buf(const double* data_ptr,
    const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count);

// Reallocating of standard-layout matrix to
// double block representation with b1 x b2 main blocks
// and d1 x d2 subblocks, using a buffer array.
double* standard_to_double_block_layout_reallocation_buf(
    const double* data_ptr,
    const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count,
    const int& d_block_rows_count, const int& d_block_cols_count);

// Inplace reallocating of standard-layout matrix to
// double block representation with b1 x b2 main blocks
// and db1 x db2 small blocks.
double* standard_to_double_block_layout_reallocation(double* data_ptr,
    const int& m_rows, const int& m_cols,
    const int& b_rows, const int& b_cols,
    const int& db_rows, const int& db_cols);

#endif  // _REALLOCATION_H_
