﻿#include "include/realloccore.h"
#include "reallocation.h"
#include "service.h"

#include <cmath>
#include <algorithm>

using std::sort;
using std::swap;
using std::min;

// *** LAYOUT DISPATCHER CLASS ***

// Controls processing of data structures, which stores such reallocation data
// as SDR-vectors. Controls memory allocation and prevents memory leakage.
// To use reallocation functions user must call 'InitDispatchSystem' function.
// To free memory, which is occupied by contained structs, user must call
// 'TurfOffDispatchSystem' function after using reallocation functions.
// Methods 'init' and 'clear' aren't intended for user utilization.
// Method 'CreateNewRecord' creates object of class 'BlockReallocationInfo'.
// Method 'find' finds struct corresponding to specified layout parameters.
class LayoutDataDispatcher
{
public:
    static vector<BlockReallocationInfo*> cache_data;
    static BlockReallocationInfo trivial_layout_info;

    static void init()
    {
        cache_data.clear();
    }

    static void clear()
    {
        for (size_t i = 0; i < cache_data.size(); ++i)
        {
            const_cast<BlockReallocationInfo*>(cache_data[i])->~BlockReallocationInfo();
        }
        cache_data.clear();
    }

    static const BlockReallocationInfo* createNewRecord(
        const vector<int>& sdr_main,
        const vector<int>& sdr_addit,
        const TaskClass& main_data,
        const TaskClass& main_data_addit)
    {
        BlockReallocationInfo* new_record = new BlockReallocationInfo();
        new_record->main_data = main_data;
        new_record->main_data_addit = main_data_addit;
        new_record->sdr_main = sdr_main;
        new_record->sdr_addit = sdr_addit;

        cache_data.push_back(new_record);


        /*printf("\n Adding %d %d %d %d \n",
            main_data.getDataRef().M_ROWS,
            main_data.getDataRef().M_COLS,
            main_data.getDataRef().B_ROWS,
            main_data.getDataRef().B_COLS);*/

        return new_record;
    }

    static const BlockReallocationInfo* find(
        const int N1, const int N2,
        const int B1, const int B2)
    {
        // Value for trivial layout
        if (N2 == 0) return &trivial_layout_info;
        if (B2 == 0) return &trivial_layout_info;
        if (N1 == 0) return &trivial_layout_info;
        if (B1 == 0) return &trivial_layout_info;

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

vector<BlockReallocationInfo*> LayoutDataDispatcher::cache_data =
                                    vector<BlockReallocationInfo*>();
BlockReallocationInfo LayoutDataDispatcher::trivial_layout_info =
                                    BlockReallocationInfo();

void InitDispatchSystem()
{
    LayoutDataDispatcher::init();
}

void TurnOffDispatchSystem()
{
    LayoutDataDispatcher::trivial_layout_info.~BlockReallocationInfo();
    LayoutDataDispatcher::clear();
}

// *** LAYOUT DISPATCHER CLASS ***


double* map_with_block_layout(
    double* dst_ptr,
    const double* src_ptr,
    const TaskClass& task_info)
{
    const TaskData& data = task_info.getDataRef();
    const int size = data.M_ROWS * data.M_COLS;
    for (int i = 0; i < size; ++i)
    {
        dst_ptr[task_info.indexFunction(i)] = src_ptr[i];
    }
    return dst_ptr;
}

double* map_with_double_block_layout(
    double* dst_ptr,
    const double* src_ptr,
    const TaskClass& task_info)
{
    const TaskData& data = task_info.getDataRef();
    for (int i = 0; i < data.M_ROWS; ++i)
    {
        const double* src_row = src_ptr + i * data.M_COLS;
        for (int j = 0; j < data.M_COLS; ++j)
        {
            dst_ptr[task_info.indexFunctionDbl(i,j)] = src_row[j];
        }
    }
    return dst_ptr;
}

// Проверка, является ли цикл, содержащий индекс 'index',
// новым по отношению к циклам, порожденным множеством индексов 'help_vec'
inline bool _fastcall is_new_cycle(const int start_index,
                                   const TaskClass& task_data,
                                   const vector<int>& help_vec)
{
    int next = start_index;
    do
    {
        // Если цикл новый, то все его индексы больше индекса,
        // с которого начинали поиск, т.к. поиск производится
        // слева направо. Поэтому, встречая индекс меньший
        // стартового, можно точно утверждать, что он принадлежит
        // уже пройденному циклу.
        if (next < start_index)
        {
            return false;
        }
        if (bin_search(next, help_vec))
        {
            return false;
        }
        next = task_data.indexFunctionReduced(next);
    }
    while (next != start_index);

    return true;
}

const vector<int> cycles_distribution_computation(const TaskClass& task_info)
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

    const int gcd_val = gcd(data.B_COLS, data.DIF_COLS);
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

            if (!bin_search(min_index, help_vec))
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
    sdr_vec_main = cycles_distribution_computation(main_task_info);
    // Searching for cycles of additional case (N1 % B1 != 0)
    const TaskData& task_data = main_task_info.getDataRef();
    if (task_data.DIF_ROWS != 0)
    {
        subtask_info = TaskClass(task_data.M_ROWS, task_data.M_COLS,
            task_data.DIF_ROWS, task_data.B_COLS);
        sdr_vec_addit = cycles_distribution_computation(subtask_info);
    }

    // Constructing new reallocation data struct
    const BlockReallocationInfo* new_realloc_info =
        LayoutDataDispatcher::createNewRecord(
            sdr_vec_main,
            sdr_vec_addit,
            main_task_info,
            subtask_info);

    return new_realloc_info;
}

void reallocate_stripe(double* stripe_data,
                       const TaskClass& task_info,
                       const vector<int>& sdr_vec,
                       const bool inverse = false)
{
    const TaskData& data = task_info.getDataRef();

    if ( sdr_vec.empty()      ||
        (stripe_data == NULL) ||
        (data.B_ROWS == 0)    ||
        (data.B_COLS == 0) )
    {
        return;
    }

    int(_fastcall TaskClass::*index_function) (int) const = inverse ?
        &TaskClass::indexFunctionReducedInverse :
        &TaskClass::indexFunctionReduced;

    double* buf_mem = (double*) calloc(2*data.B_COLS, sizeof(double));
    double* src_data_buffer = new(buf_mem)               double[data.B_COLS];
    double* dst_data_buffer = new(buf_mem + data.B_COLS) double[data.B_COLS];

    const int gcd_val = gcd(data.B_COLS, data.DIF_COLS);
    const int standard_width = sizeof(double) * ((data.DIF_COLS == 0) ?
        data.B_COLS : (((data.DIF_COLS > 1) && (gcd_val > 1)) ? gcd_val : 1));

    size_t sdr_vec_iterator = 0;
    const size_t sdr_vec_size = sdr_vec.size();
    // Reallocation within the bounds of current block stripe
    while (sdr_vec_iterator < sdr_vec_size)
    {
        int dst_index = sdr_vec[sdr_vec_iterator];
        int width = standard_width;
        if (++sdr_vec_iterator < sdr_vec_size)
        {
            const int next_el = sdr_vec[sdr_vec_iterator];
            // Non-standard width of cycle is stored as negative number.
            if (next_el < 0)
            {
                width = (-next_el) * sizeof(double);
                ++sdr_vec_iterator;
            }
        }

        memcpy(src_data_buffer, stripe_data + dst_index, width);
        // Holds place, where we start this cycle
        const int first_in_cycle = dst_index;
        // Do until cycle beginning isn't reached
        do
        {
            // Calculating of index for destination of copying
            dst_index = (task_info.*index_function)(dst_index);
            // Pointer to place, that is destination of copying
            double* dst_data_ptr = stripe_data + dst_index;

            // Saving of data, which is locating at destination location now
            memcpy(dst_data_buffer, dst_data_ptr, width);
            // Copying of data from source location to destination location
            memcpy(dst_data_ptr, src_data_buffer, width);
            // Swapping of source and destination buffer pointes.
            // Destination data, which was saved at this iteration,
            // will be used as source data at next iteration.
            swap(src_data_buffer, dst_data_buffer);
        }
        while (dst_index != first_in_cycle);
    }

    free(buf_mem);
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
    const TaskClass& main_task_data = realloc_info.main_data;
    const TaskClass& addit_task_data = realloc_info.main_data_addit;
    const vector<int>& sdr_vec = realloc_info.sdr_main;
    const vector<int>& sdr_vec_last_stripe = realloc_info.sdr_addit;
    
    const TaskData& data = main_task_data.getDataRef();
    double* stripe_data_ptr = data_ptr;
    // Every iteration corresponds to the separate stripe
    for (int it = 0; it < data.M_ROWS / data.B_ROWS; ++it)
    {
        reallocate_stripe(stripe_data_ptr, main_task_data, sdr_vec, true);
        stripe_data_ptr += data.STRIPE_SIZE;
    }

    reallocate_stripe(stripe_data_ptr,
                      addit_task_data,
                      sdr_vec_last_stripe,
                      true);

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

void standard_to_block_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2)
{
    // Searching for required reallocation data in cache
    const BlockReallocationInfo* realloc_info = 
        LayoutDataDispatcher::find( N1, N2, B1, B2);

    // If data was found (paramenters was already used before this call),
    // then produce reallocation without cycles searching
    if (realloc_info != NULL)
    {
        standard_to_block_layout_reallocation(data_ptr, *realloc_info);
    }
    else    // else, use standard algorithm and push results to cache
    {
        realloc_info = computeCyclesDistribution(N1, N2, B1, B2);
        standard_to_block_layout_reallocation(data_ptr, *realloc_info);
    }
}

double* block_to_standard_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2)
{
    // Searching reallocation data in cache
    const BlockReallocationInfo* realloc_info =
        LayoutDataDispatcher::find(N1, N2, B1, B2);

    // If data was found, then produce reallocation without cycles searching
    if (realloc_info != NULL)
    {
        block_to_standard_layout_reallocation(data_ptr, *realloc_info);
    }
    else    // else, produce cycles searching for this reallocation pattern
    {
        realloc_info = computeCyclesDistribution(N1, N2, B1, B2);
        block_to_standard_layout_reallocation(data_ptr, *realloc_info);
    }
    return data_ptr;
}

void standard_to_double_block_layout_reallocation(
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
        LayoutDataDispatcher::find(N1, N2, B1, B2);
    if (upper_level_realloc_info == NULL)
    {
        upper_level_realloc_info = computeCyclesDistribution(N1, N2, B1, B2);
    }
    const BlockReallocationInfo* main_realloc_info =
        LayoutDataDispatcher::find(B1, B2, D1, D2);
    if (main_realloc_info == NULL)
    {
        main_realloc_info = computeCyclesDistribution(B1, B2, D1, D2);
    }
    const BlockReallocationInfo* right_realloc_info = 
        LayoutDataDispatcher::find(B1, rb2, D1, rd2);
    if (right_realloc_info == NULL)
    {
        right_realloc_info = computeCyclesDistribution(B1, rb2, D1, rd2);
    }
    const BlockReallocationInfo* bottom_realloc_info =
        LayoutDataDispatcher::find(rb1, B2, rd1, D2);
    if (bottom_realloc_info == NULL)
    {
        bottom_realloc_info = computeCyclesDistribution(rb1, B2, rd1, D2);
    }
    const BlockReallocationInfo* corner_realloc_info =
        LayoutDataDispatcher::find(rb1, rb2, rd1, rd2);
    if (corner_realloc_info == NULL)
    {
        corner_realloc_info = computeCyclesDistribution(rb1, rb2, rd1, rd2);
    }

    DoubleBlockReallocationInfo new_realloc_info;
    new_realloc_info.upper_level_realloc_info = upper_level_realloc_info;
    new_realloc_info.main_realloc_info = main_realloc_info;
    new_realloc_info.right_realloc_info = right_realloc_info;
    new_realloc_info.bottom_realloc_info = bottom_realloc_info;
    new_realloc_info.corner_realloc_info = corner_realloc_info;

    standard_to_double_block_layout_reallocation(data_ptr, new_realloc_info);
}

double* double_block_to_standard_layout_reallocation(
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
        LayoutDataDispatcher::find(N1, N2, B1, B2);
    if (upper_level_realloc_info == NULL)
    {
        upper_level_realloc_info = computeCyclesDistribution(N1, N2, B1, B2);
    }
    const BlockReallocationInfo* main_realloc_info =
        LayoutDataDispatcher::find(B1, B2, D1, D2);
    if (main_realloc_info == NULL)
    {
        main_realloc_info = computeCyclesDistribution(B1, B2, D1, D2);
    }
    const BlockReallocationInfo* right_realloc_info =
        LayoutDataDispatcher::find(B1, rb2, D1, rd2);
    if (right_realloc_info == NULL)
    {
        right_realloc_info = computeCyclesDistribution(B1, rb2, D1, rd2);
    }
    const BlockReallocationInfo* bottom_realloc_info =
        LayoutDataDispatcher::find(rb1, B2, rd1, D2);
    if (bottom_realloc_info == NULL)
    {
        bottom_realloc_info = computeCyclesDistribution(rb1, B2, rd1, D2);
    }
    const BlockReallocationInfo* corner_realloc_info =
        LayoutDataDispatcher::find(rb1, rb2, rd1, rd2);
    if (corner_realloc_info == NULL)
    {
        corner_realloc_info = computeCyclesDistribution(rb1, rb2, rd1, rd2);
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
