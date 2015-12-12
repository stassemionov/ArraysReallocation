#include "reallocation.h"

#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

#include "service.h"

// #include "omp.h"

using std::vector;
using std::sort;
using std::swap;

// * Task paramenters
static int M_ROWS;             // matrix rows count
static int M_COLS;             // matrix columns count
static int B_ROWS;             // block rows count
static int B_COLS;             // block columns count
static int M_BLOCK_ROWS;       // число блоков в направлении столбца
static int M_BLOCK_COLS;       // число блоков в направлении строки
static int DIF_ROWS;           // смещение сетки по столбцу
static int DIF_COLS;           // смещение сетки по строке
static int B_SIZE;             // площаль блока
static int B_ROWS_SHIFT_SIZE;
static int B_COLS_SHIFT_SIZE;

inline void init_parameters(const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count)
{
    M_ROWS = rows_count;
    M_COLS = cols_count;
    B_ROWS = block_rows_count;
    B_COLS = block_cols_count;

    M_BLOCK_ROWS = static_cast<int>(ceil((double)M_ROWS / B_ROWS));
    M_BLOCK_COLS = static_cast<int>(ceil((double)M_COLS / B_COLS));
    DIF_COLS = M_COLS % B_COLS;
    DIF_ROWS = M_ROWS % B_ROWS;

    B_SIZE = B_COLS * B_ROWS;
    B_ROWS_SHIFT_SIZE = B_COLS*(B_ROWS - DIF_ROWS);
    B_COLS_SHIFT_SIZE = (DIF_COLS == 0) ? 0 : B_ROWS * (B_COLS - DIF_COLS);
}

int get_nearest_greater_divisor(const int& divident, const int& divisor)
{
    int dt = divident, dr = divisor;
    while (dt % dr != 0)
    {
        ++dr;
    }
    return dr;
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
    const int&& row_i = index / M_COLS;
    const int&& col_i = index % M_COLS;

    // координаты блока
    const int&& row_b = row_i / B_ROWS;
    const int&& col_b = col_i / B_COLS;
    
    // номер блока
    // const int I = row_b * M_BLOCK_COLS + col_b;
    // координаты элемента в блоке
    // int row_i_loc = row_i % b1;
    // int col_i_loc = col_i % b2;

    const int&& shift = ((col_b == M_BLOCK_COLS - 1) && (DIF_COLS != 0)) ?
        (row_i % B_ROWS) * DIF_COLS + (col_i % B_COLS) :
        (row_i % B_ROWS) * B_COLS   + (col_i % B_COLS);

    // Для нижних блоков нужно вычислить смещение
    // границы сетки от границы матрицы, если оно есть
    if ((row_b == M_BLOCK_ROWS - 1) && (DIF_ROWS != 0))
    {
        return (row_b * M_BLOCK_COLS + col_b)*B_SIZE + shift
            - row_b*B_COLS_SHIFT_SIZE - col_b*B_ROWS_SHIFT_SIZE;
    }
    else
    {
        return (row_b * M_BLOCK_COLS + col_b)*B_SIZE + shift
            - row_b*B_COLS_SHIFT_SIZE;
    }
}

bool m_find(const int& val, const vector<int>& vec)
{
    if (vec.empty())
    {
        return false;
    }
    if (val > vec[vec.size()-1])
    {
        return false;
    }

    int l = 0;
    int r = static_cast<int>(vec.size()) - 1;
    int m = 0;
    while (l < r)
    {
        m = (l + r) / 2;
        if (vec[m] < val)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    return vec[r] == val;
}

bool is_new_cycle(const int& index, vector<int>& help_vec)
{
    int next = index;
    vector<int> hv;
    int k = 1;
    size_t size = help_vec.size();
    do
    {
        if (m_find(next, help_vec))
        {
            return false;
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
    // Additional cycles representatives for new cycles searching speed-up.
    vector<int> help_vec;
    // Width of right block column (no multiplicity case)
    const int right_block_width = M_COLS % B_COLS;
    // Count of iterations to be done (sum of all cycles lengths)
    const int iteration_count = (right_block_width == 0) ?
        (B_ROWS*M_COLS - 2 * B_COLS) :
        (B_ROWS*M_COLS - B_COLS - right_block_width);

    int gcd_val = gcd(B_COLS, right_block_width);
    // Step of iterations
    // In case of multiplicity, every cycle width equals 'b2'.
    // It means, we can transfer 'b2' elements with one step)
    // Otherwise, cycle width equals GCD('b2','N2' mod 'b2').
    const int step = (right_block_width == 0) ? B_COLS :
        (((right_block_width > 1) && (gcd_val > 1)) ? gcd_val : 1);
    // Iterations starts with 'b2'-th element,
    // because first 'b2' elements don't need to be transfer.
    int i = B_COLS;
    // i = B_ROWS*M_COLS - ((right_block_width == 0) ? 
    // B_COLS : right_block_width);

    // Defines frequency of additional cycles representatives insertion.
    // Requires O(b1) memory, because we have O(b1) cycles on average.
    const int insertion_step = (B_ROWS*M_COLS) / (4 * B_ROWS);
    // This counter controls count of
    // additional cycles representatives insertion.
    int insert_counter = 0;

    int next_cycle_begining = i;
    int it = 0;
    // * Learning of current cycles distribution
    while (it < iteration_count)
    {
        // Do until cycle beginning isn't reached
        // (first element of current pass)
        int length = 0;
        int min_index = M_COLS*B_ROWS;
        int max_index = 0;
        const int first = i;
        do
        {
            // go to next position in this cycle
            i = f_ind(i);

            // Statistics collecting...
            if (min_index > i)
            {
                // Minimum index in current cycle
                min_index = i;
            }
            if (max_index < i)
            {
                // Maximum index in current cycle
                max_index = i;
            }
            // Current cycle length
            ++length;

            // Iteration count is bounded by count of elements,
            // which require to be transfer.
            it += step;

            // Insertion additional indices
            if ((++insert_counter) % insertion_step == 0)
            {
                if (!m_find(i, help_vec))
                {
                    help_vec.push_back(i);
                    sort(help_vec.begin(), help_vec.end());
                }
            }
        }
        while (i != first);

        // We don't need to collect indices of cycles of 1 element
        if (length > 1)
        {
            sdr_vec.push_back(min_index);
            if (!m_find(i, help_vec))
            {
                help_vec.push_back(i);
                sort(help_vec.begin(), help_vec.end());
            }
            // printf("\nMIN = %d \nMAX = %d\nLEN = %d\n",
            // min_index, max_index, length);
        }

        if (it < iteration_count)
        {
            // Next cycle searching
            do
            {
                next_cycle_begining += step;
            }
            while (!is_new_cycle(next_cycle_begining, help_vec));
            i = next_cycle_begining;
        }
    }
    // * End of learning

    // Recovering of task parameters
    if (virtual_block_height != -1)
    {
        init_parameters(M_ROWS, M_COLS, save_b_rows_value, B_COLS);
    }

    // printf("\n\n HELP %zd \n", help_vec.size());
    // printf("\n SDR  %zd \n\n", sdr_vec.size());

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
    const int right_block_width = M_COLS % B_COLS;

    const int gcd_val = gcd(B_COLS, right_block_width);
    const int len = sizeof(double)* ((right_block_width == 0) ? B_COLS :
        (((right_block_width > 1) && (gcd_val > 1)) ? gcd_val : 1));

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
            swap(buffer1, buffer2);
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
    const int last_stripe_h = M_ROWS % B_ROWS;
    if (last_stripe_h != 0)
    {
        sdr_vec_last_stripe = cycles_distribution_learning(last_stripe_h);
    }
    // double tt = omp_get_wtime();
    // Reallocation
    const int stripe_size = B_ROWS * M_COLS;
    double* stripe_data_ptr = data_ptr;
    // Every iteration corresponds to the separate stripe
    for (int it = 0; it < M_ROWS / B_ROWS; ++it)
    {
        reallocate_stripe(stripe_data_ptr, sdr_vec_main);
        stripe_data_ptr += stripe_size;
    }
    // tt = omp_get_wtime() - tt;
    // printf("\n REALLOCATION %f\n\n", tt);

    // Reallocation of last stripe elements,
    // if count of its rows is less than B_ROWS
    // (in this case, cycles destribution in last stripe and
    // in main part of array isn't the same).
    reallocate_stripe(stripe_data_ptr, sdr_vec_last_stripe, last_stripe_h);

    return data_ptr;
}
