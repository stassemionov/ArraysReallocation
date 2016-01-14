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

// Сделать!
// 1) Повыносить инварианты из циклов в последней функции и внутри индексной функции
// 2) Сократить требования к памяти в основном алгоритме
// 3) В последней функции написать в комментариях размеры рассматриваемой подматрицы для каждого случая
// 4) Подумать, как сократить огромные промежутки при поиске нового цикла
// 5) Сделать алгоритм двойного блочного размещения напрямую (а надо ??? : перемещений меньше, но обучение сложнее)

TaskData makeData(const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count,
    const int& double_block_rows_count,
    const int& double_block_cols_count)
{
    TaskData data;

    data.M_ROWS = rows_count;
    data.M_COLS = cols_count;
    data.B_ROWS = block_rows_count;
    data.B_COLS = block_cols_count;
    data.D_ROWS = double_block_rows_count;
    data.D_COLS = double_block_cols_count;

    data.M_BLOCK_ROWS = static_cast<int>(ceil(1.0 * data.M_ROWS / data.B_ROWS));
    data.M_BLOCK_COLS = static_cast<int>(ceil(1.0 * data.M_COLS / data.B_COLS));
    data.DIF_COLS = data.M_COLS % data.B_COLS;
    data.DIF_ROWS = data.M_ROWS % data.B_ROWS;

    return data;
}

inline int f_ind_double(const int& i_index,
                        const int& j_index,
                        const TaskData& data)
{    
    const int&& b_row = i_index / data.B_ROWS;
    const int&& b_col = j_index / data.B_COLS;
    
    // реальные размеры текущего блока
    const int& bm = (b_row == data.M_BLOCK_ROWS - 1) ?
        ((data.DIF_ROWS != 0) ? data.DIF_ROWS : data.B_ROWS) : data.B_ROWS;
    const int& bn = (b_col == data.M_BLOCK_COLS - 1) ?
        ((data.DIF_COLS != 0) ? data.DIF_COLS : data.B_COLS) : data.B_COLS;

    // смещение относительно больших блоков
    const int&& b_shift = b_row * data.B_ROWS * data.M_COLS +
                          bm * data.B_COLS * b_col;

    // реальные размеры маленьких блоков
    const int& db1 = (bm < data.D_ROWS) ? bm : data.D_ROWS;  
    const int& db2 = (bn < data.D_COLS) ? bn : data.D_COLS;

    const int&& d_row_count = static_cast<int>(ceil(1.0 * bm / db1));
    const int&& d_col_count = static_cast<int>(ceil(1.0 * bn / db2));

    // координаты элемента относительно большого блока, в котором он находится
    const int&& b_loc_i = i_index - b_row * data.B_ROWS;
    const int&& b_loc_j = j_index - b_col * data.B_COLS;

    // координаты малого блока относительно большого, в котором он находится
    const int&& d_row = b_loc_i / db1;
    const int&& d_col = b_loc_j / db2;

    const int& current_small_block_h = (d_row == d_row_count - 1) ?
        ((bm % db1 != 0) ? bm % db1 : db1) : db1;
    const int& current_small_block_w = (d_col == d_col_count - 1) ?
        ((bn % db2 != 0) ? bn % db2 : db2) : db2;

    const int&& d_shift = d_row * db1 * bn + d_col * current_small_block_h * db2;

    // смещение внутри малого блока
    const int&& d_loc_i = b_loc_i - d_row * db1;
    const int&& d_loc_j = b_loc_j - d_col * db2;
    const int&& loc_shift = d_loc_i * current_small_block_w + d_loc_j;
    
    return b_shift + d_shift + loc_shift;
}

static inline int f_ind(const int& index, const TaskData& data)
{
    const int&& index_i = index / data.M_COLS;
    const int&& index_j = index % data.M_COLS;

    const int&& block_i = index_i / data.B_ROWS;
    const int&& block_j = index_j / data.B_COLS;
        
    // Реальные размеры текущего блока
    const int& bm = (block_i == data.M_BLOCK_ROWS - 1) ?
        ((data.DIF_ROWS != 0) ? data.DIF_ROWS : data.B_ROWS) : data.B_ROWS;
    const int& bn = (block_j == data.M_BLOCK_COLS - 1) ?
        ((data.DIF_COLS != 0) ? data.DIF_COLS : data.B_COLS) : data.B_COLS;

    int&& block_shift = block_i * data.B_ROWS * data.M_COLS +
                        block_j * bm * data.B_COLS;
    int&& loc_shift = bn * (index_i - block_i * data.B_ROWS) +
                      index_j - block_j * data.B_COLS;

    return block_shift + loc_shift;
    
    /*
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
        (row_i % B_ROWS) * B_COLS + (col_i % B_COLS);

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
    }*/
}

double* standard_to_block_layout_reallocation_buf(const double* data_ptr,
    const TaskData& data)
{
    const int blocks_count_in_col = static_cast<int>(ceil(1.0 * data.M_ROWS / data.B_ROWS));
    const int blocks_count_in_row = static_cast<int>(ceil(1.0 * data.M_COLS / data.B_COLS));
    
    double* result = new double[data.M_ROWS * data.M_COLS];
    for (int i = 0; i < data.M_ROWS * data.M_COLS; ++i)
    {
        result[f_ind(i,data)] = data_ptr[i];
    }
    return result;
}


double* standard_to_double_block_layout_reallocation_buf(
    const double* data_ptr, const TaskData& data)
{
    double* result = new double[data.M_ROWS * data.M_COLS];
    for (int i = 0; i < data.M_ROWS; ++i)
    {
        for (int j = 0; j < data.M_COLS; ++j)
        {
            result[f_ind_double(i,j,data)] = data_ptr[i*data.M_COLS+j];
        }
    }
    return result;
}

bool is_new_cycle(const int& index, const TaskData& data, vector<int>& help_vec)
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
        next = f_ind(next, data);
    }
    while (next != index);

    return true;
}

vector<int> cycles_distribution_learning(const TaskData& data)
{
    // System of distinct representatives of cycles (SDR).
    vector<int> sdr_vec;
    // Additional cycles representatives for new cycles searching speed-up.
    vector<int> help_vec;
    // Count of iterations to be done (sum of all cycles lengths)
    const int iteration_count = (data.DIF_COLS == 0) ?
        (data.B_ROWS * data.M_COLS - 2 * data.B_COLS) :
        (data.B_ROWS * data.M_COLS - data.B_COLS - data.DIF_COLS);

    int gcd_val = gcd(data.B_COLS, data.DIF_COLS);
    // Step of iterations
    // In case of multiplicity, every cycle width equals 'b2'.
    // It means, we can transfer 'b2' elements with one step)
    // Otherwise, cycle width equals GCD('b2','N2' mod 'b2').
    const int step = (data.DIF_COLS == 0) ? data.B_COLS :
        (((data.DIF_COLS > 1) && (gcd_val > 1)) ? gcd_val : 1);
    // Iterations starts with 'b2'-th element,
    // because first 'b2' elements don't need to be transfer.
    int i = data.B_COLS;
    // i = data.B_ROWS * data.M_COLS - ((right_block_width == 0) ?
    // data.B_COLS : data.DIF_COLS);

    // Defines frequency of additional cycles representatives insertion.
    // Requires O(b1) memory, because we have O(b1) cycles on average.
    const int insertion_step = ((data.B_ROWS * data.M_COLS >= 4 * data.B_ROWS) ?
        ((data.B_ROWS * data.M_COLS) / (4 * data.B_ROWS)) : data.B_ROWS);
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
        int min_index = data.M_COLS * data.B_ROWS;
        int max_index = 0;
        const int first = i;
        do
        {
            // go to next position in this cycle
            i = f_ind(i, data);

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
            while (!is_new_cycle(next_cycle_begining, data, help_vec));
            i = next_cycle_begining;
        }
    }
    // * End of learning

    // printf("\n\n HELP %zd \n", help_vec.size());
    // printf("\n SDR  %zd \n\n", sdr_vec.size());

    return sdr_vec;
}

void reallocate_stripe(double* stripe_data,
                       const TaskData& data,
                       const vector<int>& sdr_vec)
{
    if (sdr_vec.empty() || (stripe_data == NULL) ||
        (data.B_ROWS == 0) || (data.B_COLS == 0))
    {
        return;
    }

    double* buffer1 = new double[data.B_COLS];
    double* buffer2 = new double[data.B_COLS];
    double* loc_data_ptr = NULL;

    const int gcd_val = gcd(data.B_COLS, data.DIF_COLS);
    const int len = sizeof(double)* ((data.DIF_COLS == 0) ? data.B_COLS :
        (((data.DIF_COLS > 1) && (gcd_val > 1)) ? gcd_val : 1));

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
            i = f_ind(i, data);
            loc_data_ptr = stripe_data + i;

            memcpy(buffer2, loc_data_ptr, len);
            memcpy(loc_data_ptr, buffer1, len);
            swap(buffer1, buffer2);
        }
        while (first_in_cycle != i);
    }
    delete[] buffer1;
    delete[] buffer2;
}

double* standard_to_block_layout_reallocation(double* data_ptr,
                                              const TaskData& data)
{
    // Systems of distinct representatives for different cycles
    vector<int> sdr_vec_main, sdr_vec_last_stripe;

    // Learning
    sdr_vec_main = cycles_distribution_learning(data);

    TaskData subtask_data;
    if (data.DIF_ROWS != 0)
    {
        subtask_data = makeData(data.M_ROWS, data.M_COLS, 
                                data.DIF_ROWS, data.B_COLS);
        sdr_vec_last_stripe = cycles_distribution_learning(subtask_data);
    }
    // double tt = omp_get_wtime();
    // Reallocation
    const int stripe_size = data.B_ROWS * data.M_COLS;
    double* stripe_data_ptr = data_ptr;
    // Every iteration corresponds to the separate stripe
    for (int it = 0; it < data.M_ROWS / data.B_ROWS; ++it)
    {
        reallocate_stripe(stripe_data_ptr, data, sdr_vec_main);
        stripe_data_ptr += stripe_size;
    }
    // tt = omp_get_wtime() - tt;
    // printf("\n REALLOCATION %f\n\n", tt);

    // Reallocation of last stripe elements,
    // if count of its rows is less than B_ROWS
    // (in this case, cycles destribution in last stripe and
    // in main part of array isn't the same).
    reallocate_stripe(stripe_data_ptr, subtask_data, sdr_vec_last_stripe);

    return data_ptr;
}


double* standard_to_double_block_layout_reallocation(double* data_ptr,
                                                     const TaskData& data)
{
    // Firstly, make block reallocation
    standard_to_block_layout_reallocation(data_ptr, data);

    // Secondly, every block must be reallocated locally
    // as well as whole matrix was reallocated just now.

    // Learning for all cases...
    vector<int> sdr_main, sdr_right, sdr_bottom, sdr_corner;
    vector<int> sdr_main_addit, sdr_right_addit,
                sdr_bottom_addit, sdr_corner_addit;
    
    // Task parameters for different learning cases.
    TaskData main_data, main_data_addit, 
             bottom_data, bottom_data_addit, 
             right_data, right_data_addit, 
             corner_data, corner_data_addit;

    // основная группа блоков (полноценные блоки)
    main_data = makeData(data.B_ROWS, data.B_COLS, data.D_ROWS, data.D_COLS);
    sdr_main = cycles_distribution_learning(main_data);
    if (data.B_ROWS % data.D_ROWS != 0)  // если локальная сетка некратна размеру основного блока
    {
        const int dif_rows_loc = (data.B_ROWS % data.D_ROWS >= data.D_ROWS) ? data.D_ROWS : data.B_ROWS % data.D_ROWS;
        main_data_addit = makeData(dif_rows_loc, data.B_COLS, dif_rows_loc, data.D_COLS);
        sdr_main_addit = cycles_distribution_learning(main_data_addit);
    }

    // группа нижних неполных блоков (усечены снизу)
    if (data.DIF_ROWS != 0)  // если глобальная сетка некратна размеру матрицы по столбцам
    {
        const int loc_block_height = (data.DIF_ROWS >= data.D_ROWS) ? data.D_ROWS : data.DIF_ROWS;
        bottom_data = makeData(data.DIF_ROWS, data.B_COLS, loc_block_height, data.D_COLS);
        sdr_bottom = cycles_distribution_learning(bottom_data);

        const int dif_rows_loc = data.DIF_ROWS % loc_block_height;
        if (dif_rows_loc != 0)  // если локальная сетка некратна размеру нижнего неполного блока
        {
            const int db1_bottom = (dif_rows_loc >= data.D_ROWS) ? data.D_ROWS : dif_rows_loc;
            bottom_data_addit = makeData(dif_rows_loc, data.B_COLS, db1_bottom, data.D_COLS);
            sdr_bottom_addit = cycles_distribution_learning(bottom_data_addit);
        }
    }

    // группа правых неполных блоков (усечены справа)
    if (data.DIF_COLS != 0) // если глобальная сетка некратна размеру матрицы по строкам
    {
        const int loc_block_width = (data.DIF_COLS >= data.D_COLS) ? data.D_COLS : data.DIF_COLS;
        right_data = makeData(data.B_ROWS, data.DIF_COLS, data.D_ROWS, loc_block_width);
        sdr_right = cycles_distribution_learning(right_data);
        const int dif_rows_loc = data.B_ROWS % data.D_ROWS;
        if (dif_rows_loc != 0)  // если локальная сетка некратна размеру правого неполного блока
        {
            const int db1_bottom = (dif_rows_loc >= data.D_ROWS) ? data.D_ROWS : dif_rows_loc;
            right_data_addit = makeData(dif_rows_loc, data.DIF_COLS, db1_bottom, loc_block_width);
            sdr_right_addit = cycles_distribution_learning(right_data_addit);
        }
    }

    // угловой блок (усечен справа и снизу)
    if ((data.DIF_ROWS != 0) && (data.DIF_COLS != 0))   // если глобальная сетка некратна размеру матрицы ни по строкам, ни по столбцам
    {
        const int loc_block_height = (data.DIF_ROWS >= data.D_ROWS) ? data.D_ROWS : data.DIF_ROWS;
        const int loc_block_width  = (data.DIF_COLS >= data.D_COLS) ? data.D_COLS : data.DIF_COLS;
        corner_data = makeData(data.DIF_ROWS, data.DIF_COLS, loc_block_height, loc_block_width);
        sdr_corner = cycles_distribution_learning(corner_data);
        const int dif_rows_loc = data.DIF_ROWS % loc_block_height;
        if (dif_rows_loc != 0)
        {
            const int db1_bottom = (dif_rows_loc >= data.D_ROWS) ? data.D_ROWS : dif_rows_loc;
            corner_data_addit = makeData(dif_rows_loc, data.DIF_COLS, db1_bottom, loc_block_width);
            sdr_corner_addit = cycles_distribution_learning(corner_data_addit);
        }
    }
    // ... end of learning.

    // Reallocation ...
    for (int ib = 0; ib < data.M_BLOCK_ROWS; ++ib)
    {
        const int b1  = ((data.DIF_ROWS != 0) && (ib == data.M_BLOCK_ROWS - 1)) ? data.DIF_ROWS : data.B_ROWS;
        const int db1 = ((data.DIF_ROWS != 0) && (ib == data.M_BLOCK_ROWS - 1) && (data.DIF_ROWS < data.D_ROWS)) ? data.DIF_ROWS : data.D_ROWS;

        for (int jb = 0; jb < data.M_BLOCK_COLS; ++jb)
        {            
            const int b2  = ((data.DIF_COLS != 0) && (jb == data.M_BLOCK_COLS - 1)) ? data.DIF_COLS : data.B_COLS;
            const int loc_stripe_size = db1 * b2;

            // Cycles distribution for main part of block, which is being reallocated
            vector<int>& sdr_vec = ((data.DIF_ROWS != 0) && (ib == data.M_BLOCK_ROWS - 1)) ?
                (((jb == data.M_BLOCK_COLS - 1) && (data.DIF_COLS != 0)) ? sdr_corner : sdr_bottom) :
                (((jb == data.M_BLOCK_COLS - 1) && (data.DIF_COLS != 0)) ? sdr_right  : sdr_main);

            // Cycles distribution for bottom part of block, which is being reallocated.
            // This part is caused by inconsistency between block size and mesh of small blocks
            vector<int>& sdr_vec_addit = ((data.DIF_ROWS != 0) && (ib == data.M_BLOCK_ROWS - 1)) ?
                (((jb == data.M_BLOCK_COLS - 1) && (data.DIF_COLS != 0)) ? sdr_corner_addit : sdr_bottom_addit) :
                (((jb == data.M_BLOCK_COLS - 1) && (data.DIF_COLS != 0)) ? sdr_right_addit  : sdr_main_addit);

            TaskData& main_parameters = ((data.DIF_ROWS != 0) && (ib == data.M_BLOCK_ROWS - 1)) ?
                (((jb == data.M_BLOCK_COLS - 1) && (data.DIF_COLS != 0)) ? corner_data : bottom_data) :
                (((jb == data.M_BLOCK_COLS - 1) && (data.DIF_COLS != 0)) ? right_data  : main_data);
            
            TaskData& addit_parameters = ((data.DIF_ROWS != 0) && (ib == data.M_BLOCK_ROWS - 1)) ?
                (((jb == data.M_BLOCK_COLS - 1) && (data.DIF_COLS != 0)) ? corner_data_addit : bottom_data_addit) :
                (((jb == data.M_BLOCK_COLS - 1) && (data.DIF_COLS != 0)) ? right_data_addit : main_data_addit);

            // указатель на начало блока
            double* stripe_data_ptr = data_ptr + ib * data.B_ROWS * data.M_COLS + jb * b1 * data.B_COLS;

            for (int it = 0; it < b1 / db1; ++it)
            {
                reallocate_stripe(stripe_data_ptr, main_parameters, sdr_vec);
                stripe_data_ptr += loc_stripe_size;
            }
            reallocate_stripe(stripe_data_ptr, addit_parameters, sdr_vec_addit);
        }
    }

    return data_ptr;
}
