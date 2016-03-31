#include "reallocation.h"
#include "service.h"

#include <cmath>
#include <algorithm>

using std::sort;
using std::swap;
using std::min;

// Сделать!
// 1) Сократить просмотр больших промежутков при поиске нового цикла
// 2) Определить более надежное условие параллельности циклов (да и это не подводит..)
// 3) Найти способ не заносить элементы параллельных циклов в help_vec
// 6) Подумать и поэкспериментировать, можно ли,
//    имея sdr для [b1,b2,d1,d2] и для [b2,b3,d2,d3], получить sdr для [b1,b3,d1,d3] ?
// 7) Применить экономичные структуры данных для оптимизации объема памяти для хранения посещенных элементов

double* standard_to_block_layout_reallocation_buf(const double* data_ptr,
    const TaskClass& task_info)
{
    const TaskData& data = task_info.getDataRef();

    const int blocks_count_in_col = static_cast<int>(ceil(1.0 * data.M_ROWS / data.B_ROWS));
    const int blocks_count_in_row = static_cast<int>(ceil(1.0 * data.M_COLS / data.B_COLS));
    
    double* result = new double[data.M_ROWS * data.M_COLS];
    for (int i = 0; i < data.M_ROWS * data.M_COLS; ++i)
    {
        result[task_info.indexFunction(i)] = data_ptr[i];
    }
    return result;
}

double* standard_to_double_block_layout_reallocation_buf(
    const double* data_ptr, const TaskClass& task_info)
{
    const TaskData& data = task_info.getDataRef();

    double* result = new double[data.M_ROWS * data.M_COLS];
    for (int i = 0; i < data.M_ROWS; ++i)
    {
        const double* row_ptr = data_ptr + i * data.M_COLS;
        for (int j = 0; j < data.M_COLS; ++j)
        {
            result[task_info.indexFunctionDbl(i,j)] = row_ptr[j];
        }
    }
    return result;
}

// Проверка, является ли цикл, содержащий индекс 'index',
// новым по отношению к циклам, порожденным множеством индексов 'help_vec'
static inline bool is_new_cycle(const int index,
                                const TaskClass& task_info,
                                const vector<int>& help_vec)
{
    int next = index;
    do
    {
        if (m_find(next, help_vec))
        {
            return false;
        }
        next = task_info.indexFunctionReduced(next);
    }
    while (next != index);

    return true;
}

const vector<int> cycles_distribution_learning(const TaskClass& task_info)
{
    const TaskData& data = task_info.getDataRef();

    // System of distinct representatives of cycles (SDR).
    // This vector has strucure:
    // <el_0> <width_0> <el_1> <width_1> ... , where
    // 'el_i' is element of SDR,
    // 'width_i' is width of cycle, which contains element 'el_i'.
    vector<int> sdr_vec;
    // Additional cycles representatives for new cycles searching speed-up.
    vector<int> help_vec;
    help_vec.reserve(data.B_ROWS);
    // Count of iterations to be done (sum of all cycles lengths)
    const int iteration_count = (data.DIF_COLS == 0) ?
        (data.B_ROWS * data.M_COLS - 2 * data.B_COLS) :
        (data.B_ROWS * data.M_COLS - data.B_COLS - data.DIF_COLS);

    int gcd_val = gcd(data.B_COLS, data.DIF_COLS);
    // Iteration step.
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
    
    int previous_min = -step - 1;
    int previous_max = -step - 1;
    int previous_len = -step - 1;
    int next_cycle_begining = i;
    int it = 0;
    // * Learning of current cycles distribution
    while (it < iteration_count)
    {
        // Do until cycle beginning isn't reached
        // (first element of current pass)
        int length = 0;
        int inserted = 0;
        int min_index = data.M_COLS * data.B_ROWS;
        int max_index = 0;
        const int first = i;
        do
        {
            // go to next position in this cycle
            i = task_info.indexFunctionReduced(i);

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
            if (++insert_counter == insertion_step)
            {
                help_vec.push_back(i);
                ++inserted;
                insert_counter = 0;
            }
        }
        while (i != first);

        // Sort for binary search.
        if (inserted > 0)
        {
            sort(help_vec.begin(), help_vec.end());
        }

        // We don't need to collect indices of cycles with 1 element
        if (length > 1)
        {
            bool is_the_same_cycle = false;
            // Если новый цикл параллелен предыдущему,
            // то это части одного широкого цикла
            if (min_index == previous_min + step)
            if (max_index == previous_max + step)
            if (length == previous_len)
            if (task_info.indexFunctionReduced(min_index) == task_info.indexFunctionReduced(previous_min) + step)
            {
                is_the_same_cycle = true;
                // Если текущий цикл уже отмечен как цикл нестандартной ширины
                if (sdr_vec.back() < 0)
                {
                    // то его ширина увеличивается на 'step'
                    sdr_vec.back() -= step;
                }
                else
                {
                    // Иначе, обнаружен цикл, параллельный текущему;
                    // обозначаем нестандартную ширину цикла -
                    // двойную от стандартной
                    sdr_vec.push_back(-2*step);
                }
            }
            if (!is_the_same_cycle)
            {
                sdr_vec.push_back(min_index);
            }
            previous_min = min_index;
            previous_max = max_index;
            previous_len = length;

            if (!m_find(min_index, help_vec))
            {
                help_vec.push_back(min_index);
                sort(help_vec.begin(), help_vec.end());
            }

  //          printf("\nMIN = %d\nMAX = %d\nLEN = %d\n",
  //          min_index, max_index, length);
            //for (int uu = 0; uu < sdr_vec.size(); ++uu)
            //    printf("%d ", sdr_vec[uu]);
            //printf("\n\n");
            //system("pause");
        }

        if (it < iteration_count)
        {
            // Next cycle searching
            do
            {
                next_cycle_begining += step;
            }
            while (!is_new_cycle(next_cycle_begining, task_info, help_vec));
            i = next_cycle_begining;
        }
    }
    // * End of learning
    
  //   printf("\n\n HELP %zd\n", help_vec.size());
  //   printf("\n SDR %zd\n\n", sdr_vec.size());
  //   system("pause");

    return sdr_vec;
}

static const BlockReallocationInfo* computeCyclesDistribution(
    const int N1, const int N2,
    const int B1, const int B2)
{
    TaskClass main_task_info(N1, N2, B1, B2);
    TaskClass subtask_info;
    vector<int> sdr_vec_main;
    vector<int> sdr_vec_addit;

    // Searching for cycles
    sdr_vec_main = cycles_distribution_learning(main_task_info);
    // Searching for cycles of additional case (N1 % B1 != 0)
    const TaskData& task_data = main_task_info.getDataRef();
    if (task_data.DIF_ROWS != 0)
    {
        subtask_info.makeData(task_data.M_ROWS, task_data.M_COLS,
            task_data.DIF_ROWS, task_data.B_COLS);
        sdr_vec_addit = cycles_distribution_learning(subtask_info);
    }

    // Constructing of new reallocation data struct
    BlockReallocationInfo* new_realloc_info = new BlockReallocationInfo;
    new_realloc_info->main_data = main_task_info;
    new_realloc_info->main_data_addit = subtask_info;
    new_realloc_info->sdr_main = sdr_vec_main;
    new_realloc_info->sdr_addit = sdr_vec_addit;

    return new_realloc_info;
}

void reallocate_stripe(double* stripe_data,
                       const TaskClass& task_info,
                       const vector<int>& sdr_vec)
{
    const TaskData& data = task_info.getDataRef();

    if (sdr_vec.empty() || (stripe_data == NULL) ||
        (data.B_ROWS == 0) || (data.B_COLS == 0))
    {
        return;
    }

    double* buffer1 = new double[data.B_COLS];
    double* buffer2 = new double[data.B_COLS];
    double* loc_data_ptr = NULL;

    const int gcd_val = gcd(data.B_COLS, data.DIF_COLS);
    const int standard_width = sizeof(double) * ((data.DIF_COLS == 0) ?
        data.B_COLS : (((data.DIF_COLS > 1) && (gcd_val > 1)) ? gcd_val : 1));

    size_t cycle_counter = 0;
    // Reallocation within the bounds of current block stripe
    while (cycle_counter < sdr_vec.size())
    {
        int i = sdr_vec[cycle_counter];
        int width = standard_width;
        if (++cycle_counter < sdr_vec.size())
        {
            int next_el = sdr_vec[cycle_counter];
            // Non-standard width of cycle is stored as negative number.
            if (next_el < 0)
            {
                width = (-next_el) * sizeof(double);
                ++cycle_counter;
            }
        }

        memcpy(buffer1, stripe_data + i, width);
        // Remember place, where we start this cycle
        const int first_in_cycle = i;
        // Do until cycle beginning isn't reached
        do
        {
            i = task_info.indexFunctionReduced(i);
            loc_data_ptr = stripe_data + i;

            memcpy(buffer2, loc_data_ptr, width);
            memcpy(loc_data_ptr, buffer1, width);
            swap(buffer1, buffer2);
        }
        while (first_in_cycle != i);
    }
    delete[] buffer1;
    delete[] buffer2;
}

void reallocate_stripe_inverse(double* stripe_data,
                               const TaskClass& task_info,
                               const vector<int>& sdr_vec)
{
    const TaskData& data = task_info.getDataRef();

    if (sdr_vec.empty() || (stripe_data == NULL) ||
        (data.B_ROWS == 0) || (data.B_COLS == 0))
    {
        return;
    }

    double* buffer1 = new double[data.B_COLS];
    double* buffer2 = new double[data.B_COLS];
    double* loc_data_ptr = NULL;

    const int gcd_val = gcd(data.B_COLS, data.DIF_COLS);
    const int standard_width = sizeof(double) * ((data.DIF_COLS == 0) ?
        data.B_COLS : (((data.DIF_COLS > 1) && (gcd_val > 1)) ? gcd_val : 1));

    size_t cycle_counter = 0;
    // Reallocation within the bounds of current block stripe
    while (cycle_counter < sdr_vec.size())
    {
        int i = sdr_vec[cycle_counter];
        int width = standard_width;
        if (++cycle_counter < sdr_vec.size())
        {
            int next_el = sdr_vec[cycle_counter];
            // Non-standard width of cycle is stored as negative number.
            if (next_el < 0)
            {
                width = (-next_el) * sizeof(double);
                ++cycle_counter;
            }
        }

        memcpy(buffer1, stripe_data + i, width);
        // Remember place, where we start this cycle
        const int first_in_cycle = i;
        // Do until cycle beginning isn't reached
        do
        {
            i = task_info.indexFunctionReducedInverse(i);
            loc_data_ptr = stripe_data + i;

            memcpy(buffer2, loc_data_ptr, width);
            memcpy(loc_data_ptr, buffer1, width);
            swap(buffer1, buffer2);
        } while (first_in_cycle != i);
    }
    delete[] buffer1;
    delete[] buffer2;
}

const BlockReallocationInfo* standard_to_block_layout_reallocation(
                                                double* data_ptr,
                                                const TaskClass& task_info)
{
    const TaskData& data = task_info.getDataRef();

    // Systems of distinct representatives for different cycles
    vector<int> sdr_vec_main, sdr_vec_last_stripe;

    // Learning
    sdr_vec_main = cycles_distribution_learning(task_info);

    TaskClass subtask_info;
    if (data.DIF_ROWS != 0)
    {
        subtask_info.makeData(data.M_ROWS, data.M_COLS, 
                              data.DIF_ROWS, data.B_COLS);
        sdr_vec_last_stripe = cycles_distribution_learning(subtask_info);
    }
    // Reallocation
    double* stripe_data_ptr = data_ptr;
    // Every iteration corresponds to the separate stripe
    for (int it = 0; it < data.M_ROWS / data.B_ROWS; ++it)
    {
        reallocate_stripe(stripe_data_ptr, task_info, sdr_vec_main);
        stripe_data_ptr += data.STRIPE_SIZE;
    }

    // Reallocation of last stripe elements,
    // if count of its rows is less than B_ROWS
    // (in this case, cycles destribution in last stripe and
    // in main part of array isn't the same).
    reallocate_stripe(stripe_data_ptr, subtask_info, sdr_vec_last_stripe);

    BlockReallocationInfo* res_info = new BlockReallocationInfo;
    res_info->main_data = task_info.getDataRef();
    res_info->main_data_addit = subtask_info;
    res_info->sdr_main = sdr_vec_main;
    res_info->sdr_addit = sdr_vec_last_stripe;

    return res_info;
}

void standard_to_block_layout_reallocation(
    double* data_ptr,
    const BlockReallocationInfo& realloc_info)
{
    const TaskClass& main_task_data = realloc_info.main_data;
    const TaskClass& subtask_data = realloc_info.main_data_addit;
    const vector<int>& sdr_vec_main = realloc_info.sdr_main;
    const vector<int>& sdr_vec_addit = realloc_info.sdr_addit;

    double* stripe_data_ptr = data_ptr;
    const TaskData& main_task_params = main_task_data.getDataRef();
    const int it_count = main_task_params.M_ROWS / main_task_params.B_ROWS;
    // Every iteration corresponds to the separate stripe
    for (int it = 0; it < it_count; ++it)
    {
        reallocate_stripe(stripe_data_ptr, main_task_data, sdr_vec_main);
        stripe_data_ptr += main_task_params.STRIPE_SIZE;
    }
    reallocate_stripe(stripe_data_ptr, subtask_data, sdr_vec_addit);
}


const DoubleBlockReallocationInfo* standard_to_double_block_layout_reallocation(
                                                     double* data_ptr,
                                                     const TaskClass& task_info)
{
    DoubleBlockReallocationInfo* realloc_info = new DoubleBlockReallocationInfo;

    // Firstly, make block reallocation
    realloc_info->upper_level_realloc_info =
        standard_to_block_layout_reallocation(data_ptr, task_info);

    // Secondly, every block must be reallocated locally
    // as well as whole matrix was reallocated just now.

    const TaskData& data = task_info.getDataRef();

    // Learning for all cases...
    vector<int> sdr_main, sdr_right, sdr_bottom, sdr_corner;
    vector<int> sdr_main_addit, sdr_right_addit,
                sdr_bottom_addit, sdr_corner_addit;
    
    // Task parameters for different learning cases.
    TaskClass main_data, main_data_addit, 
              bottom_data, bottom_data_addit, 
              right_data, right_data_addit, 
              corner_data, corner_data_addit;

    // Learning for different relationships of small and big blocks is produced below.
    // Two main cases are highlighted:
    // 1. Learning of small blocks distribution in main part of big block.
    // 2. Learning of small blocks distribution in truncated part of big block.
    // Big block is b1 x b2 block that can be truncated 
    // on right and/or bottom side, when there is no multiplicity
    // of matrix size and block size.
    // Truncated part is part of big block, which is defined by remainder of
    // division big block size by small block size.
    // Main part is part of big block, which is defined by
    // cutting truncated part off from big block.

    // Learning for main group of big blocks (consists of complete blocks)
    main_data.makeData(data.B_ROWS, data.B_COLS, data.D_ROWS, data.D_COLS);
    sdr_main = cycles_distribution_learning(main_data);
    const int loc_rows_shift = data.B_ROWS % data.D_ROWS;
    // If there is no multiplicity of big and small blocks sizes
    if (loc_rows_shift != 0)
    {
        const int dif_rows_loc = (loc_rows_shift >= data.D_ROWS) ? 
                                  data.D_ROWS : loc_rows_shift;
        main_data_addit.makeData(dif_rows_loc, data.B_COLS,
                                 dif_rows_loc, data.D_COLS);
        sdr_main_addit = cycles_distribution_learning(main_data_addit);
    }

    // Learning for bottom group of big blocks
    // (consists of blocks, which is truncated on bottom side)
    if (data.DIF_ROWS != 0)  // если глобальная сетка некратна размеру матрицы по столбцам
    {
        const int loc_block_height = (data.DIF_ROWS >= data.D_ROWS) ?
                                      data.D_ROWS : data.DIF_ROWS;
        bottom_data.makeData(data.DIF_ROWS, data.B_COLS,
                             loc_block_height, data.D_COLS);
        sdr_bottom = cycles_distribution_learning(bottom_data);

        const int dif_rows_loc = data.DIF_ROWS % loc_block_height;
        // If there is no multiplicity of bottom big and small block size
        if (dif_rows_loc != 0)
        {
            const int db1_bottom = (dif_rows_loc >= data.D_ROWS) ?
                                    data.D_ROWS : dif_rows_loc;
            bottom_data_addit.makeData(dif_rows_loc, data.B_COLS,
                                       db1_bottom, data.D_COLS);
            sdr_bottom_addit = cycles_distribution_learning(bottom_data_addit);
        }
    }
    
    // Learning for right group of big blocks
    // (consists of blocks, which is truncated on right side)
    if (data.DIF_COLS != 0) // если глобальная сетка некратна размеру матрицы по строкам
    {
        const int loc_block_width = (data.DIF_COLS >= data.D_COLS) ?
                                     data.D_COLS : data.DIF_COLS;
        right_data.makeData(data.B_ROWS, data.DIF_COLS,
                            data.D_ROWS, loc_block_width);
        sdr_right = cycles_distribution_learning(right_data);
        const int dif_rows_loc = loc_rows_shift;
        // If there is no multiplicity of right big and small block size
        if (dif_rows_loc != 0)
        {
            const int db1_bottom = (dif_rows_loc >= data.D_ROWS) ?
                                    data.D_ROWS : dif_rows_loc;
            right_data_addit.makeData(dif_rows_loc, data.DIF_COLS,
                                      db1_bottom, loc_block_width);
            sdr_right_addit = cycles_distribution_learning(right_data_addit);
        }
    }

    // Learning for corner big block 
    // (this block can be truncated on right and bottom side)
    if ((data.DIF_ROWS != 0) && (data.DIF_COLS != 0))
    {
        const int loc_block_height = (data.DIF_ROWS >= data.D_ROWS) ?
                                      data.D_ROWS : data.DIF_ROWS;
        const int loc_block_width  = (data.DIF_COLS >= data.D_COLS) ?
                                      data.D_COLS : data.DIF_COLS;
        corner_data.makeData(data.DIF_ROWS, data.DIF_COLS,
                             loc_block_height, loc_block_width);
        sdr_corner = cycles_distribution_learning(corner_data);
        const int dif_rows_loc = data.DIF_ROWS % loc_block_height;
        // If there is no multiplicity of corner big and small block size
        if (dif_rows_loc != 0)
        {
            const int db1_bottom = (dif_rows_loc >= data.D_ROWS) ? 
                                    data.D_ROWS : dif_rows_loc;
            corner_data_addit.makeData(dif_rows_loc, data.DIF_COLS,
                                       db1_bottom, loc_block_width);
            sdr_corner_addit = cycles_distribution_learning(corner_data_addit);
        }
    }
    // ... end of learning.

    BlockReallocationInfo* main_realloc_info = new BlockReallocationInfo;
    BlockReallocationInfo* right_realloc_info = new BlockReallocationInfo;
    BlockReallocationInfo* bottom_realloc_info = new BlockReallocationInfo;
    BlockReallocationInfo* corner_realloc_info = new BlockReallocationInfo;

    main_realloc_info->main_data = main_data;
    main_realloc_info->main_data_addit = main_data_addit;
    main_realloc_info->sdr_main = sdr_main;
    main_realloc_info->sdr_addit = sdr_main_addit;

    right_realloc_info->main_data = right_data;
    right_realloc_info->main_data_addit = right_data_addit;
    right_realloc_info->sdr_main = sdr_right;
    right_realloc_info->sdr_addit = sdr_right_addit;

    bottom_realloc_info->main_data = bottom_data;
    bottom_realloc_info->main_data_addit = bottom_data_addit;
    bottom_realloc_info->sdr_main = sdr_bottom;
    bottom_realloc_info->sdr_addit = sdr_bottom_addit;

    corner_realloc_info->main_data = corner_data;
    corner_realloc_info->main_data_addit = corner_data_addit;
    corner_realloc_info->sdr_main = sdr_corner;
    corner_realloc_info->sdr_addit = sdr_corner_addit;

    realloc_info->main_realloc_info = main_realloc_info;
    realloc_info->right_realloc_info = right_realloc_info;
    realloc_info->bottom_realloc_info = bottom_realloc_info;
    realloc_info->corner_realloc_info = corner_realloc_info;

    // Reallocation ...
    for (int ib = 0; ib < data.M_BLOCK_ROWS; ++ib)
    {
        const bool is_bottom_incomplete_stripe = (data.DIF_ROWS != 0) &&
                                                 (ib == data.M_BLOCK_ROWS - 1);

        const int b1  = is_bottom_incomplete_stripe ? data.DIF_ROWS : data.B_ROWS;
        const int db1 = (is_bottom_incomplete_stripe &&
                        (data.DIF_ROWS < data.D_ROWS)) ?
                         data.DIF_ROWS : data.D_ROWS;

        for (int jb = 0; jb < data.M_BLOCK_COLS; ++jb)
        {
            const bool is_right_incomplete_stripe = (data.DIF_COLS != 0) &&
                                                    (jb == data.M_BLOCK_COLS - 1);

            const int b2  = is_right_incomplete_stripe ? data.DIF_COLS :
                                                         data.B_COLS;
            const int loc_stripe_size = db1 * b2;

            // Cycles distribution for main part of block,
            // which is being reallocated
            vector<int>* sdr_vec;

            // Cycles distribution for truncated part of block,
            // which is being reallocated.
            vector<int>* sdr_vec_addit;

            TaskClass *main_parameters, *addit_parameters;
            if (is_bottom_incomplete_stripe)
            {
                if (is_right_incomplete_stripe)
                {
                    sdr_vec = &sdr_corner;
                    sdr_vec_addit = &sdr_corner_addit;
                    main_parameters = &corner_data;
                    addit_parameters = &corner_data_addit;
                }
                else
                {
                    sdr_vec = &sdr_bottom;
                    sdr_vec_addit = &sdr_bottom_addit;
                    main_parameters = &bottom_data;
                    addit_parameters = &bottom_data_addit;
                }
            }
            else
            {
                if (is_right_incomplete_stripe)
                {
                    sdr_vec = &sdr_right;
                    sdr_vec_addit = &sdr_right_addit;
                    main_parameters = &right_data;
                    addit_parameters = &right_data_addit;
                }
                else
                {
                    sdr_vec = &sdr_main;
                    sdr_vec_addit = &sdr_main_addit;
                    main_parameters = &main_data;
                    addit_parameters = &main_data_addit;
                }
            }

            // Pointer on first element of big block which is being reallocated
            double* stripe_data_ptr = data_ptr +
                                      ib * data.B_ROWS * data.M_COLS +
                                      jb * b1 * data.B_COLS;

            // Reallocation of main part of big block
            for (int it = 0; it < b1 / db1; ++it)
            {
                reallocate_stripe(stripe_data_ptr,
                                  *main_parameters,
                                  *sdr_vec);
                stripe_data_ptr += loc_stripe_size;
            }
            // Reallocation of truncated part of big block
            reallocate_stripe(stripe_data_ptr,
                              *addit_parameters,
                              *sdr_vec_addit);
        }
    }

    return realloc_info;
}

void standard_to_double_block_layout_reallocation(
    double* data_ptr,
    const DoubleBlockReallocationInfo& realloc_info)
{
    // Firstly, make block reallocation
    standard_to_block_layout_reallocation(data_ptr, *realloc_info.upper_level_realloc_info);

    // Secondly, every block must be reallocated locally
    // as well as whole matrix was reallocated just now.

    // Unpacking reallocation data from parameter 'realloc_info'
    const TaskData& nlevel_data = realloc_info.upper_level_realloc_info->main_data.getDataRef();
    const TaskData& blevel_data = realloc_info.main_realloc_info->main_data.getDataRef();
    // SDR-vectors for different cycles distribution cases
    const vector<int>& sdr_main = realloc_info.main_realloc_info->sdr_main;
    const vector<int>& sdr_main_addit = realloc_info.main_realloc_info->sdr_addit;
    const vector<int>& sdr_right = realloc_info.right_realloc_info->sdr_main;
    const vector<int>& sdr_right_addit = realloc_info.right_realloc_info->sdr_addit;
    const vector<int>& sdr_bottom = realloc_info.bottom_realloc_info->sdr_main;
    const vector<int>& sdr_bottom_addit = realloc_info.bottom_realloc_info->sdr_addit;
    const vector<int>& sdr_corner = realloc_info.corner_realloc_info->sdr_main;
    const vector<int>& sdr_corner_addit = realloc_info.corner_realloc_info->sdr_addit;
    // Task parameters for different cycles distribution cases
    const TaskClass& main_data = realloc_info.main_realloc_info->main_data;
    const TaskClass& main_data_addit = realloc_info.main_realloc_info->main_data_addit;
    const TaskClass& right_data = realloc_info.right_realloc_info->main_data;
    const TaskClass& right_data_addit = realloc_info.right_realloc_info->main_data_addit;
    const TaskClass& bottom_data = realloc_info.bottom_realloc_info->main_data;
    const TaskClass& bottom_data_addit = realloc_info.bottom_realloc_info->main_data_addit;
    const TaskClass& corner_data = realloc_info.corner_realloc_info->main_data;
    const TaskClass& corner_data_addit = realloc_info.corner_realloc_info->main_data_addit;

    // Reallocation ...
    for (int ib = 0; ib < nlevel_data.M_BLOCK_ROWS; ++ib)
    {
        const bool is_bottom_incomplete_stripe = (nlevel_data.DIF_ROWS != 0) &&
            (ib == nlevel_data.M_BLOCK_ROWS - 1);

        const int b1 = is_bottom_incomplete_stripe ? nlevel_data.DIF_ROWS : nlevel_data.B_ROWS;
        const int db1 = (is_bottom_incomplete_stripe &&
            (nlevel_data.DIF_ROWS < blevel_data.B_ROWS)) ?
            nlevel_data.DIF_ROWS : blevel_data.B_ROWS;

        for (int jb = 0; jb < nlevel_data.M_BLOCK_COLS; ++jb)
        {
            const bool is_right_incomplete_stripe = (nlevel_data.DIF_COLS != 0) &&
                (jb == nlevel_data.M_BLOCK_COLS - 1);

            const int b2 = is_right_incomplete_stripe ? nlevel_data.DIF_COLS :
                nlevel_data.B_COLS;
            const int loc_stripe_size = db1 * b2;

            // Cycles distribution for main part of block,
            // which is being reallocated
            const vector<int>* sdr_vec;
            // Cycles distribution for truncated part of block,
            // which is being reallocated.
            const vector<int>* sdr_vec_addit;
            const TaskClass *main_parameters, *addit_parameters;

            if (is_bottom_incomplete_stripe)
            {
                if (is_right_incomplete_stripe)
                {
                    sdr_vec = &sdr_corner;
                    sdr_vec_addit = &sdr_corner_addit;
                    main_parameters = &corner_data;
                    addit_parameters = &corner_data_addit;
                }
                else
                {
                    sdr_vec = &sdr_bottom;
                    sdr_vec_addit = &sdr_bottom_addit;
                    main_parameters = &bottom_data;
                    addit_parameters = &bottom_data_addit;
                }
            }
            else
            {
                if (is_right_incomplete_stripe)
                {
                    sdr_vec = &sdr_right;
                    sdr_vec_addit = &sdr_right_addit;
                    main_parameters = &right_data;
                    addit_parameters = &right_data_addit;
                }
                else
                {
                    sdr_vec = &sdr_main;
                    sdr_vec_addit = &sdr_main_addit;
                    main_parameters = &main_data;
                    addit_parameters = &main_data_addit;
                }
            }

            // Pointer on first element of big block which is being reallocated
            double* stripe_data_ptr = data_ptr +
                ib * nlevel_data.STRIPE_SIZE +
                jb * b1 * nlevel_data.B_COLS;

            // Reallocation of main part of big block
            for (int it = 0; it < b1 / db1; ++it)
            {
                reallocate_stripe(stripe_data_ptr,
                    *main_parameters,
                    *sdr_vec);
                stripe_data_ptr += loc_stripe_size;
            }
            // Reallocation of truncated part of big block
            reallocate_stripe(stripe_data_ptr,
                *addit_parameters,
                *sdr_vec_addit);
        }
    }
}

double* block_to_standard_layout_reallocation(double* data_ptr,
                                              const BlockReallocationInfo& realloc_info)
{
    const TaskClass main_task_data = realloc_info.main_data;
    const TaskClass addit_task_data = realloc_info.main_data_addit;
    const vector<int>& sdr_vec = realloc_info.sdr_main;
    const vector<int>& sdr_vec_last_stripe = realloc_info.sdr_addit;
    
    const TaskData& data = main_task_data.getDataRef();
    double* stripe_data_ptr = data_ptr;
    // Every iteration corresponds to the separate stripe
    for (int it = 0; it < data.M_ROWS / data.B_ROWS; ++it)
    {
        reallocate_stripe_inverse(stripe_data_ptr, main_task_data, sdr_vec);
        stripe_data_ptr += data.STRIPE_SIZE;
    }

    reallocate_stripe_inverse(stripe_data_ptr, addit_task_data, sdr_vec_last_stripe);

    return data_ptr;
}

double* double_block_to_standard_layout_reallocation(
                            double* data_ptr,
                            const DoubleBlockReallocationInfo& realloc_info)
{
    const BlockReallocationInfo* block_realloc_info = realloc_info.upper_level_realloc_info;
    const TaskData& block_realloc_params = block_realloc_info->main_data.getDataRef();
    const int N1 = block_realloc_params.M_ROWS;
    const int N2 = block_realloc_params.M_COLS;
    const int NB1 = block_realloc_params.M_BLOCK_ROWS;
    const int NB2 = block_realloc_params.M_BLOCK_COLS;
    const int B1 = block_realloc_params.B_ROWS;
    const int B2 = block_realloc_params.B_COLS;

    // Reallocates each big block with inverse block reallocation algorythm
    const BlockReallocationInfo* local_realloc_info = NULL;
    for (int ib = 0; ib < NB1; ++ib)
    {
        const int row_shift = ib*B1*N2;
        const int rb1 = min(B1, N1 - ib*B1);
        for (int jb = 0; jb < NB2; ++jb)
        {
            const int rb2 = min(B2, N2 - jb*B2);
            double* block_begin = data_ptr + row_shift + jb*rb1*B2;

            if (rb1 == B1)
            {
                if (rb2 == B2)
                    local_realloc_info = realloc_info.main_realloc_info;
                else
                    local_realloc_info = realloc_info.right_realloc_info;
            }
            else
            {
                if (rb2 == B2)
                    local_realloc_info = realloc_info.bottom_realloc_info;
                else
                    local_realloc_info = realloc_info.corner_realloc_info;
            }
            block_to_standard_layout_reallocation(block_begin, *local_realloc_info);
        }
    }
    // Reallocates whole matrix after each block reallocation
    block_to_standard_layout_reallocation(data_ptr, *block_realloc_info);
    return data_ptr;
}




// * RELEASE VERSIONS * //

class ReallocationDataCache
{
public:
    static vector<const BlockReallocationInfo*> cache_data;

    static void add(const BlockReallocationInfo* record)
    {
        cache_data.push_back(record);
    }

    static void clear()
    {
        cache_data.clear();
    }

    static const BlockReallocationInfo* find(
        const int N1, const int N2,
        const int B1, const int B2)
    {
        for (size_t i = 0; i < cache_data.size(); ++i)
        {
            const BlockReallocationInfo* info = cache_data[i];
            const TaskData& data = info->main_data.getDataRef();
            if (data.M_ROWS == N1)
            {
                if (data.M_COLS == N2)
                {
                    if (data.B_ROWS == B1)
                    {
                        if (data.B_COLS == B2)
                        {
                            return info;
                        }
                    }
                }
            }
        }
        return NULL;
    }
};

vector<const BlockReallocationInfo*> ReallocationDataCache::cache_data =
    vector<const BlockReallocationInfo*>();


void standard_to_block_layout_reallocation_release(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2)
{
    // Searching reallocation data in cache
    const BlockReallocationInfo* realloc_info = 
        ReallocationDataCache::find( N1, N2, B1, B2);

    // If there (paramenters was used before this call),
    // then produce reallocation without cycles searching
    if (realloc_info != NULL)
    {
        standard_to_block_layout_reallocation(data_ptr, *realloc_info);
    }
    else    // else, use standard algorithm and push results to cache
    {
        TaskClass params(N1, N2, B1, B2);
        const BlockReallocationInfo* new_realloc_info =
            standard_to_block_layout_reallocation(data_ptr, params);
        ReallocationDataCache::add(new_realloc_info);
    }
}

double* block_to_standard_layout_reallocation_release(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2)
{
    // Searching reallocation data in cache
    const BlockReallocationInfo* realloc_info =
        ReallocationDataCache::find(N1, N2, B1, B2);

    // If there, then produce inverse reallocation without cycles searching
    if (realloc_info != NULL)
    {
        block_to_standard_layout_reallocation(data_ptr, *realloc_info);
    }
    else    // else, produce cycles searching for this reallocation pattern
    {
        const BlockReallocationInfo* new_realloc_info =
            computeCyclesDistribution(N1, N2, B1, B2);
        // Inverse reallocation
        block_to_standard_layout_reallocation(data_ptr, *new_realloc_info);
        // Pushing created struct to cache
        ReallocationDataCache::add(new_realloc_info);
    }
    return data_ptr;
}

void standard_to_double_block_layout_reallocation_release(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2)
{
    const int rb1 = N1 % B1;
    const int rb2 = N2 % B2;
    const int rd1 = min(D1, rb1);
    const int rd2 = min(D2, rb2);

    // Searching reallocation data for each case in cache
    const BlockReallocationInfo* upper_level_realloc_info =
        ReallocationDataCache::find(N1, N2, B1, B2);
    if (upper_level_realloc_info == NULL)
    {
        upper_level_realloc_info = computeCyclesDistribution(N1, N2, B1, B2);
        ReallocationDataCache::add(upper_level_realloc_info);
    }
    const BlockReallocationInfo* main_realloc_info =
        ReallocationDataCache::find(B1, B2, D1, D2);
    if (main_realloc_info == NULL)
    {
        main_realloc_info = computeCyclesDistribution(B1, B2, D1, D2);
        ReallocationDataCache::add(main_realloc_info);
    }
    const BlockReallocationInfo* right_realloc_info =
        ReallocationDataCache::find(B1, rb2, D1, rd2);
    if (right_realloc_info == NULL)
    {
        right_realloc_info = computeCyclesDistribution(B1, rb2, D1, rd2);
        ReallocationDataCache::add(right_realloc_info);
    }
    const BlockReallocationInfo* bottom_realloc_info =
        ReallocationDataCache::find(rb1, B2, rd1, D2);
    if (bottom_realloc_info == NULL)
    {
        bottom_realloc_info = computeCyclesDistribution(rb1, B2, rd1, D2);
        ReallocationDataCache::add(bottom_realloc_info);
    }
    const BlockReallocationInfo* corner_realloc_info =
        ReallocationDataCache::find(rb1, rb2, rd1, rd2);
    if (corner_realloc_info == NULL)
    {
        corner_realloc_info = computeCyclesDistribution(rb1, rb2, rd1, rd2);
        ReallocationDataCache::add(corner_realloc_info);
    }

    DoubleBlockReallocationInfo new_realloc_info;
    new_realloc_info.upper_level_realloc_info = upper_level_realloc_info;
    new_realloc_info.main_realloc_info = main_realloc_info;
    new_realloc_info.right_realloc_info = right_realloc_info;
    new_realloc_info.bottom_realloc_info = bottom_realloc_info;
    new_realloc_info.corner_realloc_info = corner_realloc_info;

    standard_to_double_block_layout_reallocation(data_ptr, new_realloc_info);
}

double* double_block_to_standard_layout_reallocation_release(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2)
{
    const int rb1 = N1 % B1;
    const int rb2 = N2 % B2;
    const int rd1 = min(D1, rb1);
    const int rd2 = min(D2, rb2);

    // Searching reallocation data for each case in cache
    const BlockReallocationInfo* upper_level_realloc_info =
        ReallocationDataCache::find(N1, N2, B1, B2);
    if (upper_level_realloc_info == NULL)
    {
        upper_level_realloc_info = computeCyclesDistribution(N1, N2, B1, B2);
        ReallocationDataCache::add(upper_level_realloc_info);
    }
    const BlockReallocationInfo* main_realloc_info =
        ReallocationDataCache::find(B1, B2, D1, D2);
    if (main_realloc_info == NULL)
    {
        main_realloc_info = computeCyclesDistribution(B1, B2, D1, D2);
        ReallocationDataCache::add(main_realloc_info);
    }
    const BlockReallocationInfo* right_realloc_info =
        ReallocationDataCache::find(B1, rb2, D1, rd2);
    if (right_realloc_info == NULL)
    {
        right_realloc_info = computeCyclesDistribution(B1, rb2, D1, rd2);
        ReallocationDataCache::add(right_realloc_info);
    }
    const BlockReallocationInfo* bottom_realloc_info =
        ReallocationDataCache::find(rb1, B2, rd1, D2);
    if (bottom_realloc_info == NULL)
    {
        bottom_realloc_info = computeCyclesDistribution(rb1, B2, rd1, D2);
        ReallocationDataCache::add(bottom_realloc_info);
    }
    const BlockReallocationInfo* corner_realloc_info =
        ReallocationDataCache::find(rb1, rb2, rd1, rd2);
    if (corner_realloc_info == NULL)
    {
        corner_realloc_info = computeCyclesDistribution(rb1, rb2, rd1, rd2);
        ReallocationDataCache::add(corner_realloc_info);
    }

    DoubleBlockReallocationInfo new_realloc_info;
    new_realloc_info.upper_level_realloc_info = upper_level_realloc_info;
    new_realloc_info.main_realloc_info = main_realloc_info;
    new_realloc_info.right_realloc_info = right_realloc_info;
    new_realloc_info.bottom_realloc_info = bottom_realloc_info;
    new_realloc_info.corner_realloc_info = corner_realloc_info;

    // Inverse reallocation
    double_block_to_standard_layout_reallocation(data_ptr, new_realloc_info);

    return data_ptr;
}