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

// Сделать!
// 1) Поубирать статические переменные и функции init_parameters, и вообще все статическое
// 2) Повыносить инварианты из циклов в последней функции
// 3) Сократить требования к памяти в основном алгоритме
// 4) Переделать аргументы функции reallocate_stripe - один лишний
// 5) В последней функции написать в комментариях размеры рассматриваемой подматрицы для каждого случая
// 6) Подумать, как сократить огромные промежутки при поиске нового цикла
// 7) Сделать алгоритм двойного блочного размещения напрямую (а надо ??? : перемещений меньше, но обучение сложнее)

static inline void init_parameters(const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count)
{
    M_ROWS = rows_count;
    M_COLS = cols_count;
    B_ROWS = block_rows_count;
    B_COLS = block_cols_count;

    M_BLOCK_ROWS = static_cast<int>(ceil(1.0 * M_ROWS / B_ROWS));
    M_BLOCK_COLS = static_cast<int>(ceil(1.0 * M_COLS / B_COLS));
    DIF_COLS = M_COLS % B_COLS;
    DIF_ROWS = M_ROWS % B_ROWS;

    B_SIZE = B_COLS * B_ROWS;
    B_ROWS_SHIFT_SIZE = B_COLS*(B_ROWS - DIF_ROWS);
    B_COLS_SHIFT_SIZE = (DIF_COLS == 0) ? 0 : B_ROWS * (B_COLS - DIF_COLS);
}

inline int f_ind_double(
    const int& i_index, const int& j_index,
    const int& M,  const int& N,
    const int& B1, const int& B2,
    const int& D1, const int& D2)
{
    const int&& b_row_count = static_cast<int>(ceil(1.0 * M / B1));
    const int&& b_col_count = static_cast<int>(ceil(1.0 * N / B2));
    
    // смещение относительно больших блоков
    const int&& b_row = i_index / B1;
    const int&& b_col = j_index / B2;
    
    const int& bm = (b_row == b_row_count - 1) ?    // реальные размеры текущего блока
        ((M % B1 != 0) ? M % B1 : B1) : B1;
    const int& bn = (b_col == b_col_count - 1) ?
        ((N % B2 != 0) ? N % B2 : B2) : B2;

    const int&& b_shift = b_row * B1 * N + bm * B2 * b_col;

    // смещение относительно малых блоков внутри большого блока
    const int& db1 = (bm < D1) ? bm : D1;  // реальные размеры маленьких блоков
    const int& db2 = (bn < D2) ? bn : D2;

    const int&& d_row_count = static_cast<int>(ceil(1.0 * bm / db1));
    const int&& d_col_count = static_cast<int>(ceil(1.0 * bn / db2));

    const int&& b_loc_i = i_index - b_row * B1; // координаты элемента относительно большого блока
    const int&& b_loc_j = j_index - b_col * B2; // в котором он находится

    const int&& d_row = b_loc_i / db1;  // координаты малого блока относительно большого,
    const int&& d_col = b_loc_j / db2;  // в котором он находится

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

static inline int f_ind(const int& index)
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
    }
}

double* standard_to_block_layout_reallocation_buf(const double* data_ptr,
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

double* standard_to_double_block_layout_reallocation_buf(
    const double* data_ptr,
    const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count,
    const int& d_block_rows_count, const int& d_block_cols_count)
{
    double* result = new double[cols_count*rows_count];
    for (int i = 0; i < rows_count; ++i)
    {
        for (int j = 0; j < cols_count; ++j)
        {
            result[f_ind_double(i, j, 
                rows_count, cols_count,
                block_rows_count, block_cols_count,
                d_block_rows_count, d_block_cols_count)] = data_ptr[i*cols_count+j];
        }
    }
    return result;
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

vector<int> cycles_distribution_learning(
    const int virtual_m,  const int virtual_n,
    const int virtual_b1, const int virtual_b2)
{
    // 'virtual_<...>' may be needed, when we want to get distribution
    // for task parameters, which differs from current.
    const int save_m = M_ROWS, save_n = M_COLS;
    const int save_b1 = B_ROWS, save_b2 = B_COLS;
    init_parameters(virtual_m, virtual_n, virtual_b1, virtual_b2);

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
    const int insertion_step = ((B_ROWS*M_COLS >= 4 * B_ROWS) ?
        ((B_ROWS*M_COLS) / (4 * B_ROWS)) : B_ROWS);
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
    init_parameters(save_m, save_n, save_b1, save_b2);

    // printf("\n\n HELP %zd \n", help_vec.size());
    // printf("\n SDR  %zd \n\n", sdr_vec.size());

    return sdr_vec;
}

void reallocate_stripe(double* stripe_data, const vector<int>& sdr_vec,
    const int virtual_m,  const int virtual_n,
    const int virtual_b1, const int virtual_b2)
{
    if (sdr_vec.empty() || (stripe_data == NULL) ||
        (virtual_b1 == 0) || (virtual_b2 == 0))
    {
        return;
    }

    const int save_m = M_ROWS, save_n = M_COLS;
    const int save_b1 = B_ROWS, save_b2 = B_COLS;
    init_parameters(virtual_m, virtual_n, virtual_b1, virtual_b2);

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

    init_parameters(save_m, save_n, save_b1, save_b2);
}

double* standard_to_block_layout_reallocation(double* data_ptr,
    const int& m_rows, const int& m_cols,
    const int& b_rows, const int& b_cols)
{
    // For correct computation of index mapping function
    init_parameters(m_rows, m_cols, b_rows, b_cols);

    // Systems of distinct representatives for different cycles
    vector<int> sdr_vec_main, sdr_vec_last_stripe;

    // Learning
    sdr_vec_main = cycles_distribution_learning(M_ROWS, M_COLS, B_ROWS, B_COLS);
    if (DIF_ROWS != 0)
    {
        sdr_vec_last_stripe = cycles_distribution_learning(
            M_ROWS, M_COLS, DIF_ROWS, B_COLS);
    }
    // double tt = omp_get_wtime();
    // Reallocation
    const int stripe_size = B_ROWS * M_COLS;
    double* stripe_data_ptr = data_ptr;
    // Every iteration corresponds to the separate stripe
    for (int it = 0; it < M_ROWS / B_ROWS; ++it)
    {
        reallocate_stripe(stripe_data_ptr, sdr_vec_main,
            B_ROWS, M_COLS, B_ROWS, B_COLS);
        stripe_data_ptr += stripe_size;
    }
    // tt = omp_get_wtime() - tt;
    // printf("\n REALLOCATION %f\n\n", tt);

    // Reallocation of last stripe elements,
    // if count of its rows is less than B_ROWS
    // (in this case, cycles destribution in last stripe and
    // in main part of array isn't the same).
    reallocate_stripe(stripe_data_ptr, sdr_vec_last_stripe,
        DIF_ROWS, M_COLS, DIF_ROWS, B_COLS);

    return data_ptr;
}

double* standard_to_double_block_layout_reallocation(double* data_ptr,
    const int& m_rows, const int& m_cols,
    const int& b_rows, const int& b_cols,
    const int& db_rows, const int& db_cols)
{
    // Firstly, make block reallocation
    standard_to_block_layout_reallocation(
        data_ptr, m_rows, m_cols, b_rows, b_cols);

    // Secondly, every block must be reallocated locally
    // as well as whole matrix was reallocated just now.

    // Learning for all cases...
    vector<int> sdr_main, sdr_right, sdr_bottom, sdr_corner;
    vector<int> sdr_main_addit, sdr_right_addit,
        sdr_bottom_addit, sdr_corner_addit;

    const int dif_rows_glob = m_rows % b_rows;
    const int dif_cols_glob = m_cols % b_cols;

    // основная группа блоков (полноценные блоки)
    sdr_main = cycles_distribution_learning(b_rows, b_cols, db_rows, db_cols);
    if (b_rows % db_rows != 0)  // если локальная сетка некратна размеру основного блока
    {
        const int dif_rows_loc = (b_rows % db_rows >= db_rows) ? db_rows : b_rows % db_rows;
        sdr_main_addit = cycles_distribution_learning(
            dif_rows_loc, b_cols, dif_rows_loc, db_cols);
    }

    // группа нижних неполных блоков (усечены снизу)
    if (dif_rows_glob != 0)  // если глобальная сетка некратна размеру матрицы по столбцам
    {
        const int loc_block_height = (dif_rows_glob >= db_rows) ? db_rows : dif_rows_glob;
        sdr_bottom = cycles_distribution_learning(
            dif_rows_glob, b_cols, loc_block_height, db_cols);
        const int dif_rows_loc = dif_rows_glob % loc_block_height;
        if (dif_rows_loc != 0)  // если локальная сетка некратна размеру нижнего неполного блока
        {
            const int db1_bottom = (dif_rows_loc >= db_rows) ? db_rows : dif_rows_loc;
            sdr_bottom_addit = cycles_distribution_learning(
                dif_rows_loc, b_cols, db1_bottom, db_cols);
        }
    }

    // группа правых неполных блоков (усечены справа)
    if (dif_cols_glob != 0) // если глобальная сетка некратна размеру матрицы по строкам
    {
        const int loc_block_width = (dif_cols_glob >= db_cols) ? db_cols : dif_cols_glob;
        sdr_right = cycles_distribution_learning(
            b_rows, dif_cols_glob, db_rows, loc_block_width);
        const int dif_rows_loc = b_rows % db_rows;
        if (dif_rows_loc != 0)  // если локальная сетка некратна размеру правого неполного блока
        {
            const int db1_bottom = (dif_rows_loc >= db_rows) ? db_rows : dif_rows_loc;
            sdr_right_addit = cycles_distribution_learning(
                dif_rows_loc, dif_cols_glob, db1_bottom, loc_block_width);
        }
    }

    // угловой блок (усечен справа и снизу)
    if ((dif_rows_glob != 0) && (dif_cols_glob != 0))   // если глобальная сетка некратна размеру матрицы ни по строкам, ни по столбцам
    {
        const int loc_block_height = (dif_rows_glob >= db_rows) ? db_rows : dif_rows_glob;
        const int loc_block_width  = (dif_cols_glob >= db_cols) ? db_cols : dif_cols_glob;
        sdr_corner = cycles_distribution_learning(
            dif_rows_glob, dif_cols_glob, loc_block_height, loc_block_width);
        const int dif_rows_loc = dif_rows_glob % loc_block_height;
        if (dif_rows_loc != 0)
        {
            const int db1_bottom = (dif_rows_loc >= db_rows) ? db_rows : dif_rows_loc;
            sdr_corner_addit = cycles_distribution_learning(
                dif_rows_loc, dif_cols_glob, db1_bottom, loc_block_width);
        }
    }
    // ... end of learning.

    // Reallocation ...
    const int block_rows_count = static_cast<int>(ceil(1.0 * m_rows / b_rows));
    const int block_cols_count = static_cast<int>(ceil(1.0 * m_cols / b_cols));

    for (int ib = 0; ib < block_rows_count; ++ib)
    {
        const int b1  = ((dif_rows_glob != 0) && (ib == block_rows_count - 1)) ? dif_rows_glob : b_rows;
        const int db1 = ((dif_rows_glob != 0) && (ib == block_rows_count - 1) && (dif_rows_glob < db_rows)) ? dif_rows_glob : db_rows;
        const int db1_bottom = (b1 % db1 < db1) ? b1 % db1 : db1;

        for (int jb = 0; jb < block_cols_count; ++jb)
        {            
            const int b2  = ((dif_cols_glob != 0) && (jb == block_cols_count - 1)) ? dif_cols_glob : b_cols;
            const int db2 = ((dif_cols_glob != 0) && (jb == block_cols_count - 1) && (dif_cols_glob < db_cols)) ? dif_cols_glob : db_cols;
            const int loc_stripe_size = db1 * b2;

            // Cycles distribution for main part of block, which is being reallocated
            vector<int>& sdr_vec = ((dif_rows_glob != 0) && (ib == block_rows_count - 1)) ?
                (((jb == block_cols_count - 1) && (dif_cols_glob != 0)) ? sdr_corner : sdr_bottom) :
                (((jb == block_cols_count - 1) && (dif_cols_glob != 0)) ? sdr_right  : sdr_main);

            // Cycles distribution for bottom part of block, which is being reallocated.
            // This part is caused by inconsistency between block size and mesh of small blocks
            vector<int>& sdr_vec_addit = ((dif_rows_glob != 0) && (ib == block_rows_count - 1)) ?
                (((jb == block_cols_count - 1) && (dif_cols_glob != 0)) ? sdr_corner_addit : sdr_bottom_addit) :
                (((jb == block_cols_count - 1) && (dif_cols_glob != 0)) ? sdr_right_addit  : sdr_main_addit);

            // указатель на начало блока
            double* stripe_data_ptr = data_ptr + ib*b_rows*m_cols + jb*b1*b_cols;

            for (int it = 0; it < b1 / db1; ++it)
            {
                reallocate_stripe(stripe_data_ptr, sdr_vec, b1, b2, db1, db2);
                stripe_data_ptr += loc_stripe_size;
            }
            reallocate_stripe(stripe_data_ptr, sdr_vec_addit, b1 % db1, b2, db1_bottom, db2);
        }
    }

    return data_ptr;
}
