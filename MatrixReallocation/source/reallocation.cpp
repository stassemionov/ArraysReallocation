#include "../include/realloccore.h"
#include "../include/reallocation.h"
#include "../include/service.h"

#include <cmath>
#include <algorithm>
#include <ctime>

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


// *** BUFFERED MAPPING FUNCTIONS ***
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

// *** BUFFERED MAPPING FUNCTIONS ***

// Gives basic interface to work with
// vectors, that contains SDR for cycles
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
    SDRProcessingDispatcher() {};
    
    // Joins left- and right-passing SDR-vectors into single vector,
    // that contains SDR for all cycles of current permutation.
    inline const vector<int> getCompleteSDRVec() const
    {
        vector<int> full_sdr_vec(m_sdr_vec_left.begin(), m_sdr_vec_left.end());
        full_sdr_vec.insert(full_sdr_vec.end(), m_sdr_vec_right.begin(), m_sdr_vec_right.end());
        return full_sdr_vec;
    }

    inline bool _fastcall isNewCycle(
        const int index,
        const int left_bound,
        const int right_bound,
        const TaskClass& task_data) const
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
            next_straight  = task_data.indexFunctionReduced(next_straight);
            next_backwards = task_data.indexFunctionReducedInverse(next_backwards);
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
        printf("\n\nLEFT  (%3d) : ", m_sdr_vec_left.size());
        for (size_t i = 0; i < m_sdr_vec_left.size(); ++i)
        {
            printf("%d ", m_sdr_vec_left[i]);
        }

        printf("\nRIGHT (%3d) : ", m_sdr_vec_right.size());
        for (size_t i = 0; i < m_sdr_vec_right.size(); ++i)
        {
            printf("%d ", m_sdr_vec_right[i]);
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
};

/*
// Gives basic interface to work with
// vectors, that contains SDR for cycles
// of current reallocation pattern.
class CycleIndexDispatcher
{
public:

    CycleIndexDispatcher() {};

    CycleIndexDispatcher(const int size)
    {
        m_sdr_vec.reserve(size);
    }

    // SDR-vector getter.
    inline const vector<int>& getSDRVecRef() const
    {
        return m_sdr_vec;
    }

    // Help vector getter.
    inline const vector<int>& getHelpVecRef() const
    {
        return m_help_vec;
    }

    // Binary search if specified index in help vector.
    inline bool hasIndex(const int index) const
    {
        return bin_search(index, m_help_vec);
    }

    // Checks if cycle containing 'index' is new regarding cycles,
    // that are generated by current sdr-vector indices.
    inline bool isNewCycle(
        const int start_index,
        const TaskClass& task_data) const
    {
        int next = start_index;
        do
        {
            // If there is new cycle, then all its indices are greater then
            // starting index, because searching is produced from left to right.
            // So, if there is index, that is less then staring one, then
            // it's possible to assert, that this index belongs to
            // already passed cycle.
            if (next < start_index)
            {
                return false;
            }
            if (bin_search(next, m_help_vec) > 0)
            {
                return false;
            }
            next = task_data.indexFunctionReduced(next);
        } while (next != start_index);
    
        return true;
    }

    // Accesses elements of help vector.
    inline int operator [] (const size_t i) const
    {
        return m_help_vec[i];
    }

    // Addes element to SDR-vector defining new cycle.
    inline void specifyNewCycle(const int index)
    {
        m_sdr_vec.push_back(index);
    }

    // Insertion of new index to help vector
    // with saving of elements order.
    inline void addHelpIndex(const int index)
    {
        m_help_vec.push_back(index);
        for (vector<int>::reverse_iterator it = m_help_vec.rend();
             it != m_help_vec.rbegin() - 1;
             ++it)
        {
            vector<int>::reverse_iterator next_it = it + 1;
            if (*next_it > *it)
            {
                swap(next_it, it);
            }
            else
            {
                break;
            }
        }
    }

    // Sorting of help vector.
    inline void regularize()
    {
        std::sort(m_help_vec.begin(), m_help_vec.end());
    }

    // Removing of indices, that are less then specified index,
    // that should be intepreted as starting index in next cycle pass.
    void removeUnused(const int border_index)
    {
        bool need_to_remove = false;
        size_t remove_index = 0;
        for (size_t k = 0; k < m_help_vec.size(); ++k)
        {
            if (m_help_vec[k] >= border_index)
            {
                need_to_remove = true;
                remove_index = k;
                break;
            }
        }
        if (need_to_remove)
        {
            m_help_vec.erase(
                m_help_vec.begin(),
                m_help_vec.begin() + remove_index + 1);
        }
    }

    // Writes nonstandard width of current cycle.
    inline void writeNonstandardWidth(const int iteration_step)
    {
        // If current cycle is already marked as non-standad-width one...
        if (m_sdr_vec.back() < 0)
        {
            // ...then its width must be increased
            // with value of statdard cycle.
            m_sdr_vec.back() -= iteration_step;
        }
        else
        {
            // Else, there is a parallel cycle found.
            // Indicating its nonstandard (double) width -
            // add to SDR-vector its value.
            m_sdr_vec.push_back(-2 * iteration_step);
        }
    }

private:
    // System of distinct representatives of cycles (SDR).
    // This vector has strucure:
    // <el_0> <width_0> <el_1> <width_1> ... , where
    // 'el_i' is element of SDR,
    // 'width_i' is width of cycle, which is defined by element 'el_i'.
    // Requires O(b1) memory, because we have O(b1) cycles on average.
    vector<int> m_sdr_vec;
    // Additional cycles representatives for new cycles searching speed-up.
    // Requires O(b1) memory, because we have O(b1) cycles on average.
    vector<int> m_help_vec;
};


const vector<int> cycles_distribution_computation(const TaskClass& task_info)
{
    double time_ = clock();
    const TaskData& data = task_info.getDataRef();

    CycleIndexDispatcher dispatcher(data.B_ROWS);
    // Count of iterations to be done (sum of all cycles lengths).
    const int iterations_count = (data.DIF_COLS == 0) ?
        (data.B_ROWS * data.M_COLS - 2 * data.B_COLS) :
        (data.B_ROWS * data.M_COLS - data.B_COLS - data.DIF_COLS);
    // Greatest Common Divisor of first and last block-column widths.
    const int gcd_val = gcd(data.B_COLS, data.DIF_COLS);
    // Iteration step.
    // In case of multiplicity, minimum width for all cycles equals 'b2'.
    // It means, we can transfer at least 'b2' elements with one step.
    // Otherwise, minimum cycle width equals GCD('b2','N2' mod 'b2').
    const int step = (data.DIF_COLS == 0) ? data.B_COLS : gcd_val;
    // Iterations starts with 'b2'-th element,
    // because first 'b2' elements don't need to be transfer.
    int i = data.B_COLS;
    // Defines frequency of additional cycles representatives insertion.
    const int insertion_step = data.M_COLS / 4;
    // This counter controls count of
    // additional cycles representatives insertion.
    int insert_counter = 0;
    
    int previous_min = -step - 1;
    int previous_max = -step - 1;
    int previous_len = -step - 1;
    int next_cycle_begining = i;
    int it = 0;
    // Do until all cycles aren't passed.
    while (it < iterations_count)
    {
        // Current cycle length.
        int length = 0;
        // Count of inserted additional elements of cycles.
        int inserted = 0;
        // Minimum index in current cycle.
        int min_index = data.M_COLS * data.B_ROWS;
        // Maximum index in current cycle.
        int max_index = 0;
        // Index, that is the starting point of current pass.
        const int first = i;
        // Do until cycle beginning isn't reached
        // (first element of current pass).
        do
        {
            // Going to the next position in current cycle
            i = task_info.indexFunctionReduced(i);

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
            it += step;

            // Additional indices insertion...
            if (++insert_counter == insertion_step)
            {
                dispatcher.addHelpIndex(i);
                ++inserted;
                insert_counter = 0;
            }
        }
        while (i != first);

        // Sort for binary search.
        if (inserted > 0)
        {
            dispatcher.regularize();
        }

        // We don't need to collect indices of cycles with 1 element
        if (length > 1)
        {
            // Sign of new cycle indication.
            bool is_the_same_cycle = false;
            // If new and previous cycles are parallel,
            // then they are parts of single wider cycle.
            if (length == previous_len)
            if (min_index == previous_min + step)
            if (max_index == previous_max + step)
            if (task_info.indexFunctionReduced(min_index) ==
                task_info.indexFunctionReduced(previous_min) + step)
            if (task_info.indexFunctionReduced(max_index) ==
                task_info.indexFunctionReduced(previous_max) + step)
            {
                is_the_same_cycle = true;
                dispatcher.writeNonstandardWidth(step);
            }
            // If new cycle was found, then insert
            // its starting index to the SDR-vector.
            if (!is_the_same_cycle)
            {
                dispatcher.specifyNewCycle(min_index);
                //printf("\nMIN = %d\nMAX = %d\nLEN = %d\n",
                //    min_index - data.B_COLS,
                //    (data.M_COLS*data.B_ROWS - max_index - 1) - data.DIF_COLS,
                //    length);
                //system("pause");
            }
            previous_min = min_index;
            previous_max = max_index;
            previous_len = length;

            if (!dispatcher.hasIndex(min_index))
            {
                dispatcher.addHelpIndex(min_index);
            }
            //printf("\nMIN = %d\nMAX = %d\nLEN = %d\n",
            //    min_index, max_index, length);
            //system("pause");
        }

        if (it < iterations_count)
        {
            // Next cycle searching
            do
            {
                next_cycle_begining += step;
            }
            while (!dispatcher.isNewCycle(next_cycle_begining, task_info));
            i = next_cycle_begining;
        }
    }
    // End of cycles distribution computing.
    
    // printf("\n\n HELP %zd\n", dispatcher.getHelpVecRef().size());
    // printf("\n SDR %zd\n\n", dispatcher.getSDRVecRef().size());

    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    printf("\nSEARCHING TIME = %lf\n", time_);
//    system("pause");

    return dispatcher.getSDRVecRef();
}
*/

const vector<int> cycles_distribution_computation(const TaskClass& task_info)
{
//    double time_ = clock();

    const TaskData& data = task_info.getDataRef();

    SDRProcessingDispatcher dispatcher;
    // Count of iterations to be done (sum of all cycles lengths).
    const int iterations_count = (data.DIF_COLS == 0) ?
        (data.B_ROWS * data.M_COLS - 2 * data.B_COLS) :
        (data.B_ROWS * data.M_COLS - data.B_COLS - data.DIF_COLS);
    // Greatest Common Divisor of first and last block-column widths.
    const int gcd_val = gcd(data.B_COLS, data.DIF_COLS);
    // Iteration step.
    // In case of multiplicity, minimum width for all cycles equals 'b2'.
    // It means, we can transfer at least 'b2' elements with one step.
    // Otherwise, minimum cycle width equals GCD('b2','N2' mod 'b2').
    const int step_value = (data.DIF_COLS == 0) ? data.B_COLS : gcd_val;
    // Iterations starts with 'b2'-th element,
    // because first 'b2' elements don't need to be transfer.
    int i = data.B_COLS;

    int previous_min_left = -step_value - 1;
    int previous_max_left = -step_value - 1;
    int previous_len_left = -step_value - 1;
    int next_cycle_begining_left = data.B_COLS;
    int previous_min_right = data.STRIPE_SIZE + step_value + 1;
    int previous_max_right = data.STRIPE_SIZE + step_value + 1;
    int previous_len_right = data.STRIPE_SIZE + step_value + 1;
    int next_cycle_begining_right = (data.DIF_COLS == 0) ? 
        (data.STRIPE_SIZE - data.B_COLS - 1) :
        (data.STRIPE_SIZE - data.DIF_COLS - 1);

    int it = 0;
    int direction = 1;
    int step = step_value;

    // Do until all cycles aren't passed.
    while (it < iterations_count)
    {
        // Current cycle length.
        int length = 0;
        // Minimum index in current cycle.
        int min_index = data.M_COLS * data.B_ROWS;
        // Maximum index in current cycle.
        int max_index = 0;
        // Index, that is the starting point of current pass.
        const int first = i;
        // Do until cycle beginning isn't reached
        // (first element of current pass).
        do
        {
            // Going to the next position in current cycle
            i = task_info.indexFunctionReduced(i);

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
            if (task_info.indexFunctionReduced(min_index) ==
                task_info.indexFunctionReduced(previous_min) + step)
            if (task_info.indexFunctionReduced(max_index) ==
                task_info.indexFunctionReduced(previous_max) + step)
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
                next_cycle_begining_right,
                task_info) )
            {
                starting_index += step;
            }

            i = starting_index;
        }
    }
    // End of cycles distribution computing.

//    printf("\n SDR %d\n\n", dispatcher.getStoredCount());
//    dispatcher.printSDR();

//    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
//    printf("\nSEARCHING TIME = %lf\n", time_);
//    system("pause");

    return dispatcher.getCompleteSDRVec();
}

const BlockReallocationInfo* computeCyclesDistribution(
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
