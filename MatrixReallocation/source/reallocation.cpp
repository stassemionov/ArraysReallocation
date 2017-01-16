#include "reallocation.h"
#include "service.h"

#include <cstring>
#include <algorithm>
#include <omp.h>

#include "immintrin.h"

#include <ctime>
#include <iostream>
using std::cout;
using std::endl;

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
           // const_cast<BlockReallocationInfo*>(cache_data[i])->~BlockReallocationInfo();
            delete cache_data[i];
        }
        cache_data.clear();
    }

    static const BlockReallocationInfo* createNewRecord(
        const vector<int>& sdr_main,
        const vector<int>& sdr_addit,
        const TaskClass& main_data,
        const TaskClass& main_data_addit,
        const bool transposed)
    {
        BlockReallocationInfo* new_record = new BlockReallocationInfo();
        new_record->main_data = main_data;
        new_record->main_data_addit = main_data_addit;
        new_record->sdr_main = sdr_main;
        new_record->sdr_addit = sdr_addit;
        new_record->transposed = transposed;

        cache_data.push_back(new_record);

        return new_record;
    }

    static const BlockReallocationInfo* find(
        const int N1, const int N2,
        const int B1, const int B2,
        const bool transposed = false)
    {
        // Value for trivial layout.
        if ((N2 == 0) || (B2 == 0) || (N1 == 0) || (B1 == 0))
        {
            return &trivial_layout_info;
        }
        // Value when reallocaton isn't required.
        if ((N2 == B2) && !transposed)
        {
            return createNewRecord(
                vector<int>{}, vector<int>{},
                TaskClass{N1, N2, B1, B2}, TaskClass{},
                transposed);
        }

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
                            if (info->transposed == transposed)
                            {
                                return info;
                            }
                        }
                    }
                }
            }
        }
        return nullptr;
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

// Provides basic interface to work with
// vectors, which contains SDR for cycles
// of current reallocation pattern.
class SDRProcessingDispatcher
{
    void writeWidthToLeft(const int iteration_step)
    {
        // If current cycle is already marked as non-standad-width one...
        if (m_sdr_vec_left.back() < 0)
        {
            // ...then its width must be increased
            // with value of statdard cycle.
            m_sdr_vec_left.back() -= iteration_step;
        }
        else
        {
            // Else, there is a parallel cycle found.
            // Indicating its nonstandard (double) width -
            // add to SDR-vector its value.
            m_sdr_vec_left.push_back(-2 * iteration_step);
        }
    }

    void writeWidthToRight(const int iteration_step)
    {
        int& last_el = m_sdr_vec_right.back();
        // If current cycle is already marked as non-standad-width one...
        if (last_el < 0)
        {
            // ...then its width must be increased
            // with value of statdard cycle.
            last_el += iteration_step;
            m_sdr_vec_right[m_sdr_vec_right.size() - 2] += iteration_step;
        }
        else
        {
            // Else, there is a parallel cycle found.
            // Indicating its nonstandard (double) width -
            // add to SDR-vector its value.
            last_el += iteration_step;
            m_sdr_vec_right.push_back(2 * iteration_step);
        }
    }

public:
    SDRProcessingDispatcher(const TaskClass& task_info,
        const bool transposed) :
        m_task_info{ task_info },
        m_transposed{ transposed },
        m_main_index_function{ transposed ?
            &TaskClass::indexFunctionTransposed :
            &TaskClass::indexFunctionReduced },
        m_inverse_index_function{ transposed ?
            &TaskClass::indexFunctionTransposedInverse :
            &TaskClass::indexFunctionReducedInverse }
    {};
    
    // Joins left- and right-passing SDR-vectors into single vector,
    // that contains SDR for all cycles of current permutation.
    inline const vector<int> getCompleteSDRVec() const
    {
        vector<int> full_sdr_vec(m_sdr_vec_left.begin(), m_sdr_vec_left.end());
        full_sdr_vec.insert(full_sdr_vec.end(),
            m_sdr_vec_right.begin(), m_sdr_vec_right.end());
        return full_sdr_vec;
    }

    bool isNewCycle(
        const int index,
        const int left_bound,
        const int right_bound) const
    {
        int next_straight = index;
        int next_backwards = index;
        do
        {
            if ((next_straight  < left_bound) || (next_straight  > right_bound) ||
                (next_backwards < left_bound) || (next_backwards > right_bound))
            {
                return false;
            }
            next_straight  = (m_task_info.*m_main_index_function)(next_straight);
            next_backwards = (m_task_info.*m_inverse_index_function)(next_backwards);
        }
        while ((next_straight != index) && (next_backwards != index));

        return true;
    }

    // Addes element to SDR-vector defining new cycle.
    inline void specifyNewCycle(const int index, const int direction)
    {
        if (direction > 0)
        {
            m_sdr_vec_left.push_back(index);
        }
        else
        {
            m_sdr_vec_right.push_back(index);
        }
    }

    void writeNonstandardWidth(const int iteration_step, const int direction)
    {
        if (direction > 0)
        {
            this->writeWidthToLeft(iteration_step);
        }
        else
        {
            this->writeWidthToRight(iteration_step);
        }
    }

    size_t getStoredCount(const int key = 0) const
    {
        switch (key)
        {
            case 0:
                return m_sdr_vec_left.size() + m_sdr_vec_right.size();
            case -1:
                return m_sdr_vec_right.size();
            case 1:
                return m_sdr_vec_left.size();
        }
        return 0;
    }

    void printSDR() const
    {
        printf("\n\nLEFT  (%3d) : ", (int) m_sdr_vec_left.size());
        for (size_t i = 0; i < m_sdr_vec_left.size(); ++i)
        {
            printf("%d ", m_sdr_vec_left[i]);
        }

        printf("\nRIGHT (%3d) : ", (int)m_sdr_vec_right.size());
        for (size_t i = 0; i < m_sdr_vec_right.size(); ++i)
        {
            printf("%d ", (int) m_sdr_vec_right[i]);
        }
        printf("\n\n");
    }

private:
    // System of distinct representatives of cycles (SDR).
    // These vectors have strucure as follows:
    // <el_0> <width_0> <el_1> <width_1> ... , where
    // 'el_i' is element of SDR,
    // 'width_i' is width of cycle, which is defined by element 'el_i'.
    // Require O(b1) memory, because we have O(b1) cycles on average.
    vector<int> m_sdr_vec_left;
    vector<int> m_sdr_vec_right;
    const TaskClass& m_task_info;
    const bool m_transposed;

public:
    // Mapping functions for specified layout.
    int(TaskClass::*m_main_index_function) (int) const;
    int(TaskClass::*m_inverse_index_function) (int) const;
};

const vector<int> cycles_distribution_computation(
    const TaskClass& task_info,
    const bool transposed = false)
{
//    printf(" START\n");
//    double time_ = clock();

    const TaskData& data = task_info.getDataRef();
    SDRProcessingDispatcher dispatcher(task_info, transposed);
    int(TaskClass::*index_function) (int) const =
        dispatcher.m_main_index_function;

    // Count of iterations to be done (sum of all cycles lengths).
    const int iterations_count = transposed ? (data.M_SIZE - 2) :
        ((data.DIF_COLS == 0) ? (data.STRIPE_SIZE - 2 * data.B_COLS) :
                                (data.STRIPE_SIZE - data.B_COLS - data.DIF_COLS));
    // Greatest Common Divisor of first and last block-column widths.
    const int gcd_val = gcd(data.B_COLS, data.DIF_COLS);
    // Iteration step.
    // In case of multiplicity, minimum width for all cycles equals 'b2'.
    // It means, we can transfer at least 'b2' elements with one step.
    // Otherwise, minimum cycle width equals GCD('b2','N2' mod 'b2').
    const int step_value = transposed ? 1 :
        ((data.DIF_COLS == 0) ? data.B_COLS : gcd_val);
    // Iterations starts with 'b2'-th element,
    // because first 'b2' elements don't need to be transfer.
    int i = transposed ? 1 : data.B_COLS;

    int previous_min_left = -step_value - 1;
    int previous_max_left = -step_value - 1;
    int previous_len_left = -step_value - 1;
    int next_cycle_begining_left = transposed ? 1 : data.B_COLS;
    int previous_min_right = (transposed ? data.M_SIZE : data.STRIPE_SIZE)
                            + step_value + 1;
    int previous_max_right = (transposed ? data.M_SIZE : data.STRIPE_SIZE)
                            + step_value + 1;
    int previous_len_right = (transposed ? data.M_SIZE : data.STRIPE_SIZE)
                            + step_value + 1;
    int next_cycle_begining_right = transposed ? (data.M_SIZE - 2) :
        ((data.DIF_COLS == 0) ? 
        (data.STRIPE_SIZE - data.B_COLS - 1) :
        (data.STRIPE_SIZE - data.DIF_COLS - 1));

    int it = 0;
    int direction = 1;
    int step = step_value;

    // Do until all cycles aren't passed.
    while (it < iterations_count)
    {
        // Current cycle length.
        int length = 0;
        // Minimum index in current cycle.
        int min_index = transposed ? data.M_SIZE : data.STRIPE_SIZE;
        // Maximum index in current cycle.
        int max_index = 0;
        // Index, that is the starting point of current pass.
        const int first = i;
        // Do until cycle beginning isn't reached
        // (first element of current pass).
        do
        {
            // Going to the next position in current cycle
            i = (task_info.*index_function)(i);

//            cout << i+1 << " " << iterations_count - it << endl;
//            system("pause");

            // Statistics collecting...
            if (min_index > i)
            {
                min_index = i;
            }
            if (max_index < i)
            {
                max_index = i;
            }
            ++length;
            it += step_value;
        }
        while (i != first);

        // We don't need to collect indices of cycles with 1 element
        if (length > 1)
        {
            int& previous_len = (direction > 0) ? previous_len_left : previous_len_right;
            int& previous_min = (direction > 0) ? previous_min_left : previous_min_right;
            int& previous_max = (direction > 0) ? previous_max_left : previous_max_right;
            // Sign of new cycle indication.
            bool is_the_same_cycle = false;
            // If new and previous cycles are parallel,
            // then they are parts of single wider cycle.
            if (length == previous_len)
            if (min_index == previous_min + step)
            if (max_index == previous_max + step)
            if ((task_info.*index_function)(min_index) ==
                (task_info.*index_function)(previous_min) + step)
            if ((task_info.*index_function)(max_index) ==
                (task_info.*index_function)(previous_max) + step)
            {
                is_the_same_cycle = true;
                dispatcher.writeNonstandardWidth(step, direction);
            }
            // If new cycle was found, then insert
            // its starting index to the SDR-vector.
            if (!is_the_same_cycle)
            {
                const int left_bound_of_stripe = (direction > 0) ?
                    min_index : (min_index + direction*(step_value-1));
                dispatcher.specifyNewCycle(left_bound_of_stripe, direction);
            }
            previous_len = length;
            previous_min = min_index;
            previous_max = max_index;

 //           printf(" LENGTH %d\n", length);
        }

        ((direction > 0) ?
            next_cycle_begining_left :
            next_cycle_begining_right) += step;

        if (it < iterations_count)
        {
            // Next cycle searching...
            direction = -direction;
            step = direction * step_value;
            int& starting_index = (direction > 0) ?
                next_cycle_begining_left :
                next_cycle_begining_right;

            while (!dispatcher.isNewCycle(
                starting_index,
                next_cycle_begining_left,
                next_cycle_begining_right) )
            {
                starting_index += step;
            }

            i = starting_index;
        }
    }
    // End of cycles distribution computing.

//    printf("\n SDR %d\n\n", dispatcher.getStoredCount());
//    printf(" MIN LEFT %d\n", previous_min_left);
//    printf(" MIN RIGHT %d\n",previous_min_right);
//    dispatcher.printSDR();

//    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
//    printf("\nSEARCHING TIME = %lf\n", time_);
//    system("pause");

    return dispatcher.getCompleteSDRVec();
}

// Computes generating elements for the cycles of the specified layout
// and pushs it to the cache.
const BlockReallocationInfo* computeCyclesDistribution(
    const int N1, const int N2,
    const int B1, const int B2,
    const bool transposed = false)
{
    TaskClass main_task_info(N1, N2, B1, B2);
    TaskClass subtask_info;
    vector<int> sdr_vec_main;
    vector<int> sdr_vec_addit;

    // Searching for cycles
    sdr_vec_main = cycles_distribution_computation(main_task_info, transposed);
    // Searching for cycles of additional case (N1 % B1 != 0)
    const TaskData& task_data = main_task_info.getDataRef();
    if ((task_data.DIF_ROWS != 0) && (!transposed))
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
            subtask_info,
            transposed);
    return new_realloc_info;
}

inline void copyAVX(double* dst_ptr, const double* src_ptr, const int length)
{
    int main_part = length >> 2;    // length / 4
    for (int i = 0; i < main_part; ++i)
    {
        const int shift = i << 2;   // 4 * i
        __m256d data = _mm256_set_pd(
            src_ptr[shift + 3],
            src_ptr[shift + 2],
            src_ptr[shift + 1],
            src_ptr[shift]);
    //    __m256d data = _mm256_loadu_pd(src_ptr + shift);
        _mm256_storeu_pd(dst_ptr + shift, data);
    }
    const int rest = length - ((length >> 2) << 2);
    if (rest != 0)
    {
        memcpy(dst_ptr + length - rest,
               src_ptr + length - rest,
               rest << 3);
    }
//    int main_part = length >> 3;    // length / 8
//    for (int i = 0; i < main_part; ++i)
//    {
//        const int shift = i << 3;   // 8 * i
//        __m512d data = _mm512_set_pd(
//            src_ptr[shift + 7],
//            src_ptr[shift + 6],
//            src_ptr[shift + 5],
//            src_ptr[shift + 4],
//            src_ptr[shift + 3],
//            src_ptr[shift + 2],
//            src_ptr[shift + 1],
//            src_ptr[shift] );
//        _mm512_storeu_pd(dst_ptr + shift, data);
//    }
//    const int rest = length - ((length >> 3) << 3);
//    if (rest != 0)
//    {
//        memcpy(dst_ptr + length - rest,
//            src_ptr + length - rest,
//            rest << 3);
//    }
}

// Reallocates given block stripe of matrix
// to layout specified by TaskClass parameter.
// Cycles to follow are specified by given sdr-vector.
// Direction of cycles passing is specified by 'inverse' flag.
void reallocate_stripe(double* stripe_data,
                       const TaskClass& task_info,
                       const vector<int>& sdr_vec,
                       const bool inverse = false)
{
    const TaskData& data = task_info.getDataRef();

    if (sdr_vec.empty() ||
        (stripe_data == nullptr) ||
        (data.BLOCK_SIZE == 0))
    {
        return;
    }

    int(TaskClass::*index_function) (int) const = inverse ?
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

        int w = width >> 3;
        if (w != 1)
        {
            copyAVX(src_data_buffer, stripe_data + dst_index, w);
            //        memcpy(src_data_buffer, stripe_data + dst_index, width);
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
                copyAVX(dst_data_buffer, dst_data_ptr, w);
                //            memcpy(dst_data_buffer, dst_data_ptr, width);
                // Copying of data from source location to destination location
                copyAVX(dst_data_ptr, src_data_buffer, w);
                //            memcpy(dst_data_ptr, src_data_buffer, width);
                // Swapping of source and destination buffer pointes.
                // Destination data, which was saved at this iteration,
                // will be used as source data at next iteration.
                swap(src_data_buffer, dst_data_buffer);
            }
            while (dst_index != first_in_cycle);
        }
        else
        {
            double buffer = stripe_data[dst_index];
            const int first_in_cycle = dst_index;
            do
            {
                dst_index = (task_info.*index_function)(dst_index);
                swap(buffer, stripe_data[dst_index]);
            }
            while (dst_index != first_in_cycle);
        }
    }
    free(buf_mem);
}

void reallocate_stripe_transpose(double* stripe_data,
    const TaskClass& task_info,
    const vector<int>& sdr_vec,
    const bool inverse = false)
{
    int(TaskClass::*index_function) (int) const = inverse ?
        &TaskClass::indexFunctionTransposedInverse :
        &TaskClass::indexFunctionTransposed;

    // Reallocation within the bounds of current block stripe
    for (size_t i = 0; i < sdr_vec.size(); ++i)
    {
        double src_data_buffer, dst_data_buffer;
        int dst_index = sdr_vec[i];
        src_data_buffer = stripe_data[dst_index];
        // Holds place, where we start this cycle
        const int first_in_cycle = dst_index;
        // Do until cycle beginning isn't reached
        do
        {
            // Calculating of index for destination of copying
            dst_index = (task_info.*index_function)(dst_index);
            // Saving of data, which is locating at destination location now
            dst_data_buffer = stripe_data[dst_index];
            // Copying of data from source location to destination location
            stripe_data[dst_index] = src_data_buffer;
            // Swapping of source and destination buffer pointes.
            // Destination data, which was saved at this iteration,
            // will be used as source data at next iteration.
            swap(src_data_buffer, dst_data_buffer);
        }
        while (dst_index != first_in_cycle);
    }
}

// Transposes square N x N matrix.
double* transpose_square(double* data_ptr, const int N)
{
    for (int i = 0; i < N; ++i)
    {
        double* i_row = data_ptr + i*N;
        for (int j = 0; j < i; ++j)
        {
            swap(i_row[j], data_ptr[j*N + i]);
        }
    }
    return data_ptr;
}

double* transfer_between_standard_and_block_layouts(
    double* data_ptr,
    const BlockReallocationInfo& realloc_info,
    const bool direction)
{
    const TaskClass& main_task_data = realloc_info.main_data;
    const TaskClass& subtask_data = realloc_info.main_data_addit;
    const TaskData& main_task_params = main_task_data.getDataRef();
    const vector<int>& sdr_vec_main = realloc_info.sdr_main;
    const vector<int>& sdr_vec_addit = realloc_info.sdr_addit;

    if (realloc_info.isTrivial() || (data_ptr == nullptr))
    {
        return data_ptr;
    }

    if (realloc_info.transposed)
    {
        // Condition of square transposition
        if ((main_task_params.M_ROWS == main_task_params.M_COLS) &&
            (main_task_params.M_ROWS == main_task_params.B_ROWS) &&
            (main_task_params.B_COLS == 1))
        {
            transpose_square(data_ptr, main_task_params.M_ROWS);
        }
        else
        {
            reallocate_stripe_transpose(data_ptr,
                realloc_info.main_data,
                realloc_info.sdr_main,
                !direction);
        }
        return data_ptr;
    }

    const int it_count = main_task_params.M_ROWS / main_task_params.B_ROWS;
    // Every iteration corresponds to the separate stripe.
    #pragma omp parallel for shared(main_task_data, sdr_vec_main) schedule(dynamic, 1)
    for (int it = 0; it < it_count; ++it)
    {
        reallocate_stripe(data_ptr + it * main_task_params.STRIPE_SIZE,
            main_task_data,
            sdr_vec_main,
            !direction);
    }
    reallocate_stripe(data_ptr + it_count * main_task_params.STRIPE_SIZE,
        subtask_data, sdr_vec_addit, !direction);
    return data_ptr;
}

double* standard_to_block_layout_reallocation(
    double* data_ptr,
    const BlockReallocationInfo& realloc_info)
{
    return transfer_between_standard_and_block_layouts(
        data_ptr, realloc_info, true);
}

double* block_to_standard_layout_reallocation(
    double* data_ptr,
    const BlockReallocationInfo& realloc_info)
{
    return transfer_between_standard_and_block_layouts(
        data_ptr, realloc_info, false);
}

double* transfer_between_standard_and_double_block_layouts(
    double* data_ptr,
    const DoubleBlockReallocationInfo& realloc_info,
    const bool direction)
{
    const BlockReallocationInfo* uplevel_realloc_info =
        realloc_info.upper_level_realloc_info;
    // Firstly, make block reallocation
    if (direction)
    {
        standard_to_block_layout_reallocation(
            data_ptr, *uplevel_realloc_info);
    }
    const TaskData& uplevel_data =
        uplevel_realloc_info->main_data.getDataRef();
    const int nrows = uplevel_data.M_ROWS;
    const int ncols = uplevel_data.M_COLS;
    const int brows = uplevel_data.B_ROWS;
    const int bcols = uplevel_data.B_COLS;

    // Secondly, every block must be reallocated locally
    // as well as whole matrix was reallocated just now.
    #pragma omp parallel for schedule(dynamic, 1)
    for (int ib = 0; ib < uplevel_data.M_BLOCK_ROWS; ++ib)
    {
        const int rb1 = min(brows, nrows - ib*brows);
        for (int jb = 0; jb < uplevel_data.M_BLOCK_COLS; ++jb)
        {
            const int rb2 = min(bcols, ncols - jb*bcols);
            // Cycles distribution for transposition of block,
            // which is being reallocated
            // (unused if transposition option is disabled).
            const BlockReallocationInfo* loc_realloc_info = nullptr;
            if (rb1 == brows)
            {
                loc_realloc_info = (rb2 == bcols) ?
                    realloc_info.main_realloc_info :
                    realloc_info.right_realloc_info;
            }
            else
            {
                loc_realloc_info = (rb2 == bcols) ?
                    realloc_info.bottom_realloc_info :
                    realloc_info.corner_realloc_info;
            }

            // Pointer on first element of big block which is being reallocated.
            double* block_beginning_ptr = data_ptr +
                ib * uplevel_data.STRIPE_SIZE + jb * rb1 * bcols;

            transfer_between_standard_and_block_layouts(
                block_beginning_ptr, *loc_realloc_info, direction);
        }
    }
    // Reallocates whole matrix after each block reallocation
    return direction ? data_ptr :
        transfer_between_standard_and_block_layouts(
            data_ptr, *uplevel_realloc_info, direction);
}

double* standard_to_double_block_layout_reallocation(
    double* data_ptr,
    const DoubleBlockReallocationInfo& realloc_info)
{
    return transfer_between_standard_and_double_block_layouts(
        data_ptr, realloc_info, true);
}

double* double_block_to_standard_layout_reallocation(
                            double* data_ptr,
                            const DoubleBlockReallocationInfo& realloc_info)
{
    return transfer_between_standard_and_double_block_layouts(
        data_ptr, realloc_info, false);
}


// Transpose each block in matrix with block layout,
// specified by parameters N1, N2, B1, B2.
double* transpose_each_block(double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const bool inverse = false)
{
    const int rb1 = N1 % B1;
    const int rb2 = N2 % B2;
    // Searching reallocation data for transposition.
    const BlockReallocationInfo* main_transpose_info =
        (B1 == B2) ? nullptr :
        (inverse ? LayoutDataDispatcher::find(B2, B1, B2, 1) :
                   LayoutDataDispatcher::find(B1, B2, B1, 1));
    if ((main_transpose_info == nullptr) && (B1 != B2))
    {
        main_transpose_info = inverse ?
            computeCyclesDistribution(B2, B1, B2, 1) :
            computeCyclesDistribution(B1, B2, B1, 1);
    }
    const BlockReallocationInfo* right_transpose_info =
        (B1 == rb2) ? nullptr :
        (inverse ? LayoutDataDispatcher::find(rb2, B1, rb2, 1) :
                   LayoutDataDispatcher::find(B1, rb2, B1, 1));
    if ((right_transpose_info == nullptr) && (B1 != rb2))
    {
        right_transpose_info = inverse ?
            computeCyclesDistribution(rb2, B1, rb2, 1) :
            computeCyclesDistribution(B1, rb2, B1, 1);
    }
    const BlockReallocationInfo* bottom_transpose_info =
        (rb1 == B2) ? nullptr :
        (inverse ? LayoutDataDispatcher::find(B2, rb1, B2, 1) :
                   LayoutDataDispatcher::find(rb1, B2, rb1, 1));
    if ((bottom_transpose_info == nullptr) && (rb1 != B2))
    {
        bottom_transpose_info = inverse ?
            computeCyclesDistribution(B2, rb1, B2, 1) :
            computeCyclesDistribution(rb1, B2, rb1, 1);
    }
    const BlockReallocationInfo* corner_transpose_info =
        (rb1 == rb2) ? nullptr :
        (inverse ? LayoutDataDispatcher::find(rb2, rb1, rb2, 1) :
                   LayoutDataDispatcher::find(rb1, rb2, rb1, 1));
    if ((corner_transpose_info == nullptr) && (rb1 != rb2))
    {
        corner_transpose_info = inverse ?
            computeCyclesDistribution(rb2, rb1, rb2, 1) :
            computeCyclesDistribution(rb1, rb2, rb1, 1);
    }

    const int i_ub = (int)ceil((double)N1 / B1);
    const int j_ub = (int)ceil((double)N2 / B2);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < i_ub; ++i)
    {
        const int s1 = min(B1, N1 - i*B1);
        for (int j = 0; j < j_ub; ++j)
        {
            const int s2 = min(B2, N2 - j*B2);
            double* block_beginning_ptr =
                data_ptr + i*B1*N2 + j*s1*B2;

            const BlockReallocationInfo* transpose_info = nullptr;
            if (s1 != rb1)
            {
                transpose_info = (s2 != rb2) ?
                    main_transpose_info : right_transpose_info;
            }
            else
            {
                transpose_info = (s2 != rb2) ?
                    bottom_transpose_info : corner_transpose_info;
            }

            if (transpose_info == nullptr)
            {
                transpose_square(block_beginning_ptr, s1);
            }
            else
            {
                reallocate_stripe(block_beginning_ptr,
                                  transpose_info->main_data,
                                  transpose_info->sdr_main);
            }
        }
    }
    return data_ptr;
}

// Transpose each small-block-stripe in each big block in matrix
// with block layout, specified by parameters N1, N2, B1, B2, D1.
double* transpose_stripes_in_each_block(double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1,
    const bool inverse = false)
{
    const int rb1 = N1 % B1;
    const int rb2 = N2 % B2;
    const int rd1 = B1 % D1;
    const int rbd1 = (rb1 >= D1) ? rb1 % D1 : rb1;
    // Searching reallocation data for transposition.
    const BlockReallocationInfo* main_transpose_info =
        (D1 == B2) ? nullptr :
        (inverse ? LayoutDataDispatcher::find(B2, D1, B2, 1) :
                   LayoutDataDispatcher::find(D1, B2, D1, 1));
    if ((main_transpose_info == nullptr) && (D1 != B2))
    {
        main_transpose_info = inverse ?
            computeCyclesDistribution(B2, D1, B2, 1) :
            computeCyclesDistribution(D1, B2, D1, 1);
    }
    const BlockReallocationInfo* main_addit_transpose_info =
        (rd1 == B2) ? nullptr :
        (inverse ? LayoutDataDispatcher::find(B2, rd1, B2, 1) :
                   LayoutDataDispatcher::find(rd1, B2, rd1, 1));
    if ((main_addit_transpose_info == nullptr) && (rd1 != B2))
    {
        main_addit_transpose_info = inverse ?
            computeCyclesDistribution(B2, rd1, B2, 1) :
            computeCyclesDistribution(rd1, B2, rd1, 1);
    }
    const BlockReallocationInfo* right_transpose_info =
        (D1 == rb2) ? nullptr :
        (inverse ? LayoutDataDispatcher::find(rb2, D1, rb2, 1) :
                   LayoutDataDispatcher::find(D1, rb2, D1, 1));
    if ((right_transpose_info == nullptr) && (D1 != rb2))
    {
        right_transpose_info = inverse ?
            computeCyclesDistribution(rb2, D1, rb2, 1) :
            computeCyclesDistribution(D1, rb2, D1, 1);
    }
    const BlockReallocationInfo* right_addit_transpose_info =
        (rd1 == rb2) ? nullptr :
        (inverse ? LayoutDataDispatcher::find(rb2, rd1, rb2, 1) :
                   LayoutDataDispatcher::find(rd1, rb2, rd1, 1));
    if ((right_addit_transpose_info == nullptr) && (rd1 != rb2))
    {
        right_addit_transpose_info = inverse ?
            computeCyclesDistribution(rb2, rd1, rb2, 1) :
            computeCyclesDistribution(rd1, rb2, rd1, 1);
    }
    const BlockReallocationInfo* bottom_addit_transpose_info =
        (rbd1 == B2) ? nullptr :
        (inverse ? LayoutDataDispatcher::find(B2, rbd1, B2, 1) :
                   LayoutDataDispatcher::find(rbd1, B2, rbd1, 1));
    if ((bottom_addit_transpose_info == nullptr) && (rbd1 != B2))
    {
        bottom_addit_transpose_info = inverse ?
            computeCyclesDistribution(B2, rbd1, B2, 1) :
            computeCyclesDistribution(rbd1, B2, rbd1, 1);
    }
    const BlockReallocationInfo* corner_addit_transpose_info =
        (rbd1 == rb2) ? nullptr :
        (inverse ? LayoutDataDispatcher::find(rb2, rbd1, rb2, 1) :
                   LayoutDataDispatcher::find(rbd1, rb2, rbd1, 1));
    if ((corner_addit_transpose_info == nullptr) && (rbd1 != rb2))
    {
        corner_addit_transpose_info = inverse ?
            computeCyclesDistribution(rb2, rbd1, rb2, 1) :
            computeCyclesDistribution(rbd1, rb2, rbd1, 1);
    }

    const int i_ub = (int)ceil((double)N1 / B1);  // includes last block
    const int j_ub = (int)ceil((double)N2 / B2);  // includes last block
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < i_ub; ++i)
    {
        const int s1 = min(B1, N1 - i*B1);  // rows of current block
        const int k_ub = s1 / D1;   // doesn't include last block
        const int rd = s1 - k_ub * D1;  // size of last block
        for (int j = 0; j < j_ub; ++j)
        {
            const int s2 = min(B2, N2 - j*B2);  // columns of current block
            double* block_beginning_ptr =
                data_ptr + i*B1*N2 + j*s1*B2;
            // Selection of required reallocation data for transposition.
            const BlockReallocationInfo* transpose_info = nullptr;
            const BlockReallocationInfo* transpose_info_addit = nullptr;
            if (s1 != rb1)
            {
                if (s2 != rb2)
                {
                    transpose_info = main_transpose_info;
                    transpose_info_addit = main_addit_transpose_info;
                }
                else
                {
                    transpose_info = right_transpose_info;
                    transpose_info_addit = right_addit_transpose_info;
                }
            }
            else
            {
                if (s2 != rb2)
                {
                    transpose_info = main_transpose_info;
                    transpose_info_addit = bottom_addit_transpose_info;
                }
                else
                {
                    transpose_info = right_transpose_info;
                    transpose_info_addit = corner_addit_transpose_info;
                }
            }
            // Transposition.
            if (transpose_info == nullptr)  // square small blocks
            {
                for (int k = 0; k < k_ub; ++k)
                {
                    double* stripe_ptr = block_beginning_ptr + k*D1*s2;
                    transpose_square(stripe_ptr, D1);
                }
                if (rd != 0)
                {
                    if (transpose_info_addit == nullptr)
                    {
                        transpose_square(block_beginning_ptr + k_ub*D1*s2, rd);
                    }
                    else
                    {
                        reallocate_stripe(block_beginning_ptr + k_ub*D1*s2,
                            transpose_info_addit->main_data,
                            transpose_info_addit->sdr_main);
                    }
               }
            }
            else    // nonsquare small blocks
            {
                for (int k = 0; k < k_ub; ++k)
                {
                    reallocate_stripe(block_beginning_ptr + k*D1*s2,
                        transpose_info->main_data,
                        transpose_info->sdr_main);
                }
                if (rd != 0)
                {
                    if (transpose_info_addit == nullptr)
                    {
                        transpose_square(block_beginning_ptr + k_ub*D1*s2, rd);
                    }
                    else
                    {
                        reallocate_stripe(block_beginning_ptr + k_ub*D1*s2,
                            transpose_info_addit->main_data,
                            transpose_info_addit->sdr_main);
                    }
                }
            }
        }
    }
    return data_ptr;
}




double* transfer_between_standard_and_block_layouts(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const bool direction)
{
    // If block has the same count of columns as matrix,
    // then there is no need to produce some reallocations.
    if (B2 == N2)
    {
        return data_ptr;
    }
    // Searching for required reallocation data in cache.
    const BlockReallocationInfo* realloc_info =
        LayoutDataDispatcher::find(N1, N2, B1, B2);
    // If data wasn't found (parameters wasn't used before this call),
    // then produce cycles searching and push results to cache.
    if (realloc_info == nullptr)
    {
        realloc_info = computeCyclesDistribution(N1, N2, B1, B2);
    }
    return transfer_between_standard_and_block_layouts(
        data_ptr, *realloc_info, direction);
}

double* transfer_between_standard_and_double_block_layouts(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2,
    const bool direction)
{
    // If size of main blocks equals to size of double blocks,
    // then block reallocation is enough to produce,
    // because layout of this case has trivial mapping of main block elements.
    if (B2 == D2)
    {
        return direction ?
            standard_to_block_layout_reallocation(data_ptr, N1, N2, B1, B2) :
            block_to_standard_layout_reallocation(data_ptr, N1, N2, B1, B2);
    }

    const int rb1 = N1 % B1;
    const int rb2 = N2 % B2;
    const int rd1 = min(D1, rb1);
    const int rd2 = min(D2, rb2);

    // Searching reallocation data for each case in cache
    const BlockReallocationInfo* upper_level_realloc_info =
        LayoutDataDispatcher::find(N1, N2, B1, B2);
    if (upper_level_realloc_info == nullptr)
    {
        upper_level_realloc_info = computeCyclesDistribution(N1, N2, B1, B2);
    }
    const BlockReallocationInfo* main_realloc_info =
        LayoutDataDispatcher::find(B1, B2, D1, D2);
    if (main_realloc_info == nullptr)
    {
        main_realloc_info = computeCyclesDistribution(B1, B2, D1, D2);
    }
    const BlockReallocationInfo* right_realloc_info =
        LayoutDataDispatcher::find(B1, rb2, D1, rd2);
    if (right_realloc_info == nullptr)
    {
        right_realloc_info = computeCyclesDistribution(B1, rb2, D1, rd2);
    }
    const BlockReallocationInfo* bottom_realloc_info =
        LayoutDataDispatcher::find(rb1, B2, rd1, D2);
    if (bottom_realloc_info == nullptr)
    {
        bottom_realloc_info = computeCyclesDistribution(rb1, B2, rd1, D2);
    }
    const BlockReallocationInfo* corner_realloc_info =
        LayoutDataDispatcher::find(rb1, rb2, rd1, rd2);
    if (corner_realloc_info == nullptr)
    {
        corner_realloc_info = computeCyclesDistribution(rb1, rb2, rd1, rd2);
    }

    DoubleBlockReallocationInfo new_realloc_info;
    new_realloc_info.upper_level_realloc_info = upper_level_realloc_info;
    new_realloc_info.main_realloc_info = main_realloc_info;
    new_realloc_info.right_realloc_info = right_realloc_info;
    new_realloc_info.bottom_realloc_info = bottom_realloc_info;
    new_realloc_info.corner_realloc_info = corner_realloc_info;

    return transfer_between_standard_and_double_block_layouts(
        data_ptr, new_realloc_info, direction);
}

double* transfer_between_standard_and_transposed_double_block_layouts(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2,
    const bool direction)
{
    // In this case, required layout can be obtained just with
    // usual reallocation and transposition of each big block.
    if (B1 == D1)
    {
        if (direction)
        {
            standard_to_block_layout_reallocation(data_ptr, N1, N2, B1, B2);
            transpose_each_block(data_ptr, N1, N2, B1, B2);
        }
        else
        {
            transpose_each_block(data_ptr, N1, N2, B1, B2, true);
            block_to_standard_layout_reallocation(data_ptr, N1, N2, B1, B2);
        }
        return data_ptr;
    }
    // In this case, required layout can be obtained just with
    // usual reallocation and transposition of each small block.
    if (B2 == D2)
    {
        if (direction)
        {
            standard_to_block_layout_reallocation(data_ptr, N1, N2, B1, B2);
            transpose_stripes_in_each_block(data_ptr, N1, N2, B1, B2, D1);
        }
        else
        {
            transpose_stripes_in_each_block(data_ptr, N1, N2, B1, B2, D1, true);
            block_to_standard_layout_reallocation(data_ptr, N1, N2, B1, B2);
        }
        return data_ptr;
    }

    const int rb1 = N1 % B1;
    const int rb2 = N2 % B2;
    const int rd1 = min(D1, rb1);
    const int rd2 = min(D2, rb2);

    const bool is_main_transpose = (B1 == B2) && (D1 == B1) && (D2 == 1);
    const bool is_right_transpose = (B1 == rb2) && (D1 == B1) && (rd2 == 1);
    const bool is_bottom_transpose = (rb1 == B2) && (rb1 == rd1) && (D2 == 1);
    const bool is_corner_transpose = (rb1 == rb2) && (rb1 == rd1) && (rd2 == 1);
    // Searching reallocation data for each case in cache.
    const BlockReallocationInfo* upper_level_realloc_info =
        LayoutDataDispatcher::find(N1, N2, B1, B2);
    if (upper_level_realloc_info == nullptr)
    {
        upper_level_realloc_info = computeCyclesDistribution(N1, N2, B1, B2);
    }
    const BlockReallocationInfo* main_realloc_info = is_main_transpose ?
        &LayoutDataDispatcher::trivial_layout_info :
        LayoutDataDispatcher::find(B1, B2, D1, D2, true);
    if (main_realloc_info == nullptr)
    {
        main_realloc_info = computeCyclesDistribution(B1, B2, D1, D2, true);
    }
    const BlockReallocationInfo* right_realloc_info = is_right_transpose ?
        &LayoutDataDispatcher::trivial_layout_info :
        LayoutDataDispatcher::find(B1, rb2, D1, rd2, true);
    if (right_realloc_info == nullptr)
    {
        right_realloc_info = computeCyclesDistribution(B1, rb2, D1, rd2, true);
    }
    const BlockReallocationInfo* bottom_realloc_info = is_bottom_transpose ?
        &LayoutDataDispatcher::trivial_layout_info :
        LayoutDataDispatcher::find(rb1, B2, rd1, D2, true);
    if (bottom_realloc_info == nullptr)
    {
        bottom_realloc_info = computeCyclesDistribution(rb1, B2, rd1, D2, true);
    }
    const BlockReallocationInfo* corner_realloc_info = is_corner_transpose ?
        &LayoutDataDispatcher::trivial_layout_info :
        LayoutDataDispatcher::find(rb1, rb2, rd1, rd2, true);
    if (corner_realloc_info == nullptr)
    {
        corner_realloc_info = computeCyclesDistribution(rb1, rb2, rd1, rd2, true);
    }
    DoubleBlockReallocationInfo new_realloc_info;
    new_realloc_info.upper_level_realloc_info = upper_level_realloc_info;
    new_realloc_info.main_realloc_info = main_realloc_info;
    new_realloc_info.right_realloc_info = right_realloc_info;
    new_realloc_info.bottom_realloc_info = bottom_realloc_info;
    new_realloc_info.corner_realloc_info = corner_realloc_info;

    return transfer_between_standard_and_double_block_layouts(
            data_ptr, new_realloc_info, direction);
}



// * RELEASE VERSIONS * //

double* standard_to_block_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2)
{
    return transfer_between_standard_and_block_layouts(
        data_ptr,
        N1, N2, B1, B2,
        true);
}

double* block_to_standard_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2)
{
    return transfer_between_standard_and_block_layouts(
        data_ptr,
        N1, N2, B1, B2,
        false);
}

double* standard_to_double_block_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2)
{
    return transfer_between_standard_and_double_block_layouts(
        data_ptr,
        N1, N2, B1, B2, D1, D2,
        true);
}

double* double_block_to_standard_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2)
{
    return transfer_between_standard_and_double_block_layouts(
        data_ptr,
        N1, N2, B1, B2, D1, D2,
        false);
}

double* standard_to_transposed_double_block_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2)
{
    return transfer_between_standard_and_transposed_double_block_layouts(
        data_ptr,
        N1, N2, B1, B2, D1, D2,
        true);
}

double* transposed_double_block_to_standard_layout_reallocation(
    double* data_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2)
{
    return transfer_between_standard_and_transposed_double_block_layouts(
        data_ptr,
        N1, N2, B1, B2, D1, D2,
        false);
}