#include "mappingfuncs.h"

#include <cstring>
#include <algorithm>

using std::min;

double* map_with_block_layout(
    double* dst_ptr,
    const double* src_ptr,
    const int N1, const int N2,
    const int B1, const int B2)
{
    const TaskClass task_info{ N1, N2, B1, B2 };
    for (int i = 0; i < N1 * N2; ++i)
    {
        dst_ptr[task_info.indexFunction(i)] = src_ptr[i];
    }
    return dst_ptr;
}

double* map_with_double_block_layout(
    double* dst_ptr,
    const double* src_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2)
{
    const TaskClass task_info{ N1, N2, B1, B2, D1, D2 };
    for (int i = 0; i < N1; ++i)
    {
        const double* src_row = src_ptr + i * N2;
        for (int j = 0; j < N2; ++j)
        {
            dst_ptr[task_info.indexFunctionDbl(i, j)] = src_row[j];
        }
    }
    return dst_ptr;
}

double* map_with_transposed_double_block_layout(
    double* dst_ptr,
    const double* src_ptr,
    const int N1, const int N2,
    const int B1, const int B2,
    const int D1, const int D2)
{
    const TaskClass task_info{ N1, N2, B1, B2, D1, D2 };
    for (int i = 0; i < N1 * N2; ++i)
    {
        dst_ptr[task_info.indexFunction(i)] = src_ptr[i];
    }

    int diff1 = (N1 % B1 == 0) ? B1 : N1 % B1;
    int diff2 = (N2 % B2 == 0) ? B2 : N2 % B2;
    const TaskClass main_data{ B1, B2, D1, D2 };
    const TaskClass right_data{ B1, diff2, D1, min(D2, diff2) };
    const TaskClass bottom_data{ diff1, B2, min(D1, diff1), D2 };
    const TaskClass corner_data{ diff1, diff2,
                                 min(D1, diff1), min(D2, diff2) };

    double* buffer = (double*)calloc(B1 * B2, sizeof(double));
    for (int ib = 0; ib < task_info.getDataRef().M_BLOCK_ROWS; ++ib)
    {
        const int rb1 = min(B1, N1 - ib * B1);
        for (int jb = 0; jb < task_info.getDataRef().M_BLOCK_COLS; ++jb)
        {
            const int rb2 = min(B2, N2 - jb * B2);
            const TaskClass* local_data = nullptr;
            if (rb1 == B1)
            {
                local_data = (rb2 == B2) ?
                    &main_data : &right_data;
            }
            else
            {
                local_data = (rb2 == B2) ?
                    &bottom_data : &corner_data;
            }

            // Pointer on first element of big block which is being reallocated.
            double* block_beginning_ptr = dst_ptr +
                ib * task_info.getDataRef().STRIPE_SIZE +
                jb * rb1 * B2;

            const int bsize = rb1 * rb2;

            memcpy(buffer, block_beginning_ptr, bsize * sizeof(double));

            for (int i = 0; i < bsize; ++i)
            {
                block_beginning_ptr[local_data->indexFunctionTransposed(i)] =
                    buffer[i];
            }
        }
    }
    free(buffer);
    return dst_ptr;
}