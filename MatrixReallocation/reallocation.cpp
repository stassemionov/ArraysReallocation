#include "reallocation.h"

#include <cstdlib>
#include <cmath>
#include <vector>

using std::vector;

// Task paramenters
int M_ROWS, M_COLS, B_ROWS, B_COLS;

inline void init_parameters(const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count)
{
    M_ROWS = rows_count;
    M_COLS = cols_count;
    B_ROWS = block_rows_count;
    B_COLS = block_cols_count;
}

inline int f_ind(const int& i);

double* get_reallocated(const double* data_ptr,
    const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count)
{
    // To avoid redudant argument passing to mapping function 'f_ind'
    init_parameters(rows_count, cols_count, block_rows_count, block_cols_count);

    double* result = new double[cols_count*rows_count];
    for (int i = 0; i < cols_count*rows_count; ++i)
    {
        result[f_ind(i)] = data_ptr[i];
    }
    return result;
}

inline int f_ind(const int& index)
{
    int&& row_i = static_cast<int>(index / M_COLS);
    int&& col_i = index % M_COLS;
    // координаты блока
    int&& row_b = row_i / B_ROWS;
    int&& col_b = col_i / B_COLS;

    // число блоков по столбцу (в направлении столбца)
    int&& M_b = static_cast<int>(ceil((double) M_ROWS / B_ROWS));
    // число блоков по строке (в направлении строки)
    int&& N_b = static_cast<int>(ceil((double) M_COLS / B_COLS));
    // номер блока
    int&& I = row_b * N_b + col_b;

    // координаты элемента в блоке
//  int row_i_loc = row_i % b1;
//  int col_i_loc = col_i % b2;

    int w = ( (col_b == N_b - 1) && (M_COLS % B_COLS != 0) ) ?
        M_COLS % B_COLS : B_COLS;
    int&& shift = (row_i % B_ROWS) * w + (col_i % B_COLS);
    int&& dif = M_COLS % B_COLS;

    // Для нижних блоков нужно вычислить смещение
    // границы сетки от границы матрицы
    if ((row_b == M_b - 1) && (M_ROWS % B_ROWS != 0))
    {
        if (dif == 0)
            return  I*B_COLS*B_ROWS + shift -
            col_b*B_COLS*(B_ROWS - M_ROWS % B_ROWS);
        else
            return  I*B_COLS*B_ROWS + shift -
            row_b*B_ROWS*(B_COLS - dif) -
            col_b*B_COLS*(B_ROWS - M_ROWS % B_ROWS);
    }
    else
    {
        if (dif == 0)
            return I*B_COLS*B_ROWS + shift;
        else
            return I*B_COLS*B_ROWS + shift
            - row_b*B_ROWS*(B_COLS - dif);
    }
}

bool is_new_cycle(const int& index, const vector<int>& sdr_vec)
{
    int next = index;
    size_t j, size = sdr_vec.size();
    do
    {
        for (j = 0; j < size; ++j)
        {
            if (next == sdr_vec[j])
            {
                return false;
            }
        }
        next = f_ind(next);
    }
    while (next != index);
    return true;
}

vector<int> cycles_distribution_learning(const int virtual_block_height = -1)
{
    // 'virtual_block_height' may be needed, when we want to get distribution
    // for block height, which differs from current.
    // For example, when it is necessary to get distribution in last stripe
    // when 'M_ROWS' isn't divisible by 'B_ROWS'.
    const int save_b_rows_value = B_ROWS;
    if (virtual_block_height != -1)
    {
        init_parameters(M_ROWS, M_COLS, virtual_block_height, B_COLS);
    }

    // System of distinct representatives of cycles (SDR).
    vector<int> sdr_vec;
    const int&& right_block_width = M_COLS % B_COLS;
    const int&& iteration_count = (right_block_width == 0) ?
        B_ROWS*M_COLS - 2*B_COLS : B_ROWS*M_COLS - B_COLS - right_block_width;
    int it = 0;
    const int step = (right_block_width == 0) ? B_COLS : 1;
    // Iterations starts with 'b2'-th element,
    // because first 'b2' elements don't need to be transfer.
    int i = B_COLS;
    // Maximum index among indexes of current cycle elements.
    int max_index = i;

    // * Learning of current cycles distribution
    while (it < iteration_count)
    {
        const int first_in_cycle = i;
        sdr_vec.push_back(first_in_cycle);
        // Do until cycle beginning isn't reached
        // (first element in current pass)
        do
        {
            i = f_ind(i);
            if (max_index < i)
            {
                max_index = i;
            }
            it += step;
        }
        while (first_in_cycle != i);

        if (it < iteration_count)
        {
            int&& next_cycle_begining = max_index - step;
            // Next cycle searching
            while (!is_new_cycle(next_cycle_begining, sdr_vec))
            {
                next_cycle_begining -= step;
            }
            i = next_cycle_begining;
            max_index = i;
        }
    }
    // * End of learning

    // Recovering of task parameters
    if (virtual_block_height != -1)
    {
        init_parameters(M_ROWS, M_COLS, save_b_rows_value, B_COLS);
    }

    return sdr_vec;
}

void reallocate_stripe(double* stripe_data, const vector<int>& sdr_vec,
    int virtual_block_height = -1)
{
    if ( sdr_vec.empty() || (stripe_data == NULL) )
    {
        return;
    }

    // 'virtual_block_height' may be needed, when indices in 'sdr_vec'
    // correspond distribution, which differs from current.
    const int save_b_rows_value = B_ROWS;
    if (virtual_block_height != -1)
    {
        init_parameters(M_ROWS, M_COLS, virtual_block_height, B_COLS);
    }

    double* buffer1 = new double[B_COLS];
    double* buffer2 = new double[B_COLS];
    double* loc_data_ptr;
    const int len = (M_COLS % B_COLS == 0) ?
        B_COLS * sizeof(double) : sizeof(double);

    size_t cycle_counter = 0;
    // Reallocation within the bounds of current block stripe
    while (cycle_counter < sdr_vec.size())
    {
        int i = sdr_vec[cycle_counter++];
        memcpy(buffer1, stripe_data + i, len);
        // Remember place, where we start this cycle
        const int first_in_cycle = i;
        // Do until cycle beginning isn't reached
        do
        {
            i = f_ind(i);
            loc_data_ptr = stripe_data + i;

            memcpy(buffer2, loc_data_ptr, len);
            memcpy(loc_data_ptr, buffer1, len);
            memcpy(buffer1, buffer2, len);
        }
        while (first_in_cycle != i);
    }
    delete[] buffer1;
    delete[] buffer2;

    if (virtual_block_height != -1)
    {
        init_parameters(M_ROWS, M_COLS, save_b_rows_value, B_COLS);
    }
}

double* block_reallocate_matrix(double* data_ptr, const int& m_rows,
    const int& m_cols, const int& b_rows, const int& b_cols)
{
    // For correct computation of index mapping function
    init_parameters(m_rows, m_cols, b_rows, b_cols);

    // Systems of distinct representatives for different cycles
    vector<int> sdr_vec_main, sdr_vec_last_stripe;

    // Learning
    sdr_vec_main = cycles_distribution_learning();
    const int&& last_stripe_h = M_ROWS % B_ROWS;
    if (last_stripe_h != 0)
    {
        sdr_vec_last_stripe = cycles_distribution_learning(last_stripe_h);
    }
    
    // Reallocation
    const int&& stripe_size = B_ROWS * M_COLS;
    double* stripe_data_ptr = data_ptr;
    // Every iteration corresponds to the separate stripe
    for (int it = 0; it < M_ROWS / B_ROWS; ++it)
    {
        reallocate_stripe(stripe_data_ptr, sdr_vec_main);
        stripe_data_ptr += stripe_size;
    }

    // Reallocation of last stripe elements,
    // if count of its rows is less than B_ROWS
    // (in this case, cycles destribution in last stripe and
    // in main part of array isn't the same).
    reallocate_stripe(stripe_data_ptr, sdr_vec_last_stripe, last_stripe_h);

    return data_ptr;
}
