#ifndef _REALLOCATION_H_
#define _REALLOCATION_H_

struct TaskData
{
    // Main data
    int M_ROWS;     // matrix rows count
    int M_COLS;     // matrix columns count
    int B_ROWS;     // block rows count
    int B_COLS;     // block columns count
    int D_ROWS;     // double block rows count
    int D_COLS;     // double block columns count

    // Additional data (derived from main data)
    int M_BLOCK_ROWS;       // число блоков в направлении столбца
    int M_BLOCK_COLS;       // число блоков в направлении строки
    int DIF_ROWS;           // смещение сетки по столбцу
    int DIF_COLS;           // смещение сетки по строке
};

TaskData makeData(const int& rows_count, const int& cols_count,
                  const int& block_rows_count, const int& block_cols_count,
                  const int& double_block_rows_count = 0,
                  const int& double_block_cols_count = 0);

// Inplace reallocating of standard-layout matrix to
// block representation with b1 x b2 blocks.
double* standard_to_block_layout_reallocation(double* data_ptr,
                                              const TaskData& data);

// Reallocating of standard-layout matrix to
// block representation with b1 x b2 blocks using a buffer array.
double* standard_to_block_layout_reallocation_buf(const double* data_ptr,
                                                  const TaskData& data);

// Reallocating of standard-layout matrix to
// double block representation with b1 x b2 main blocks
// and d1 x d2 subblocks, using a buffer array.
double* standard_to_double_block_layout_reallocation_buf(
    const double* data_ptr, const TaskData& data);

// Inplace reallocating of standard-layout matrix to
// double block representation with b1 x b2 main blocks
// and db1 x db2 small blocks.
double* standard_to_double_block_layout_reallocation(double* data_ptr,
                                                     const TaskData& data);

#endif  // _REALLOCATION_H_
