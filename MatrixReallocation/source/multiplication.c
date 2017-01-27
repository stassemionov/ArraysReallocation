#include <multiplication.h>

#include <stdlib.h>
#include <math.h>
#include <omp.h>

#ifdef __cplusplus
extern "C" {
#endif

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

void matrix_multiplication_double_block(double* gen_matrix,
    const double* left_matrix, const double* right_matrix,
    const int N1, const int N2, const int N3,
    const int B1, const int B2, const int B3,
    const int D1, const int D2, const int D3)
{
    const int NB1 = (int) ceil(1.0 * N1 / B1);
    const int NB2 = (int) ceil(1.0 * N2 / B2);
    const int NB3 = (int) ceil(1.0 * N3 / B3);
    const int ND1 = (int) ceil(1.0 * B1 / D1);
    const int ND2 = (int) ceil(1.0 * B2 / D2);
    const int ND3 = (int) ceil(1.0 * B3 / D3);

    // big block level
#pragma omp parallel for schedule(dynamic, 1)
    for (int ib = 0; ib < NB1; ++ib)
    {
        const int i_shift_gen = ib*B1*N3;
        const int i_shift_left = ib*B1*N2;
        const int iblock_size = min(B1, N1 - ib*B1);
        const int id_ub = min(ND1, (int) ceil(1.0 * iblock_size / D1));

        for (int jb = 0; jb < NB3; ++jb)
        {
            const int jblock_size = min(B3, N3 - jb*B3);
            const int jd_ub = min(ND3, (int) ceil(1.0 * jblock_size / D3));
            double* gen_block = gen_matrix + i_shift_gen + jb*iblock_size*B3;

            for (int kb = 0; kb < NB2; ++kb)
            {
                const int kblock_size = min(B2, N2 - kb*B2);
                const int kd_ub = min(ND2, (int) ceil(1.0 * kblock_size / D2));
                const double* left_block = left_matrix + i_shift_left + kb*iblock_size*B2;
                const double* right_block = right_matrix + kb*B2*N3 + jb*kblock_size*B3;

                // small block level
                for (int id = 0; id < id_ub; ++id)
                {
                    const int i_loc_shift = id * D1;
                    const int id_shift_gen = i_loc_shift * jblock_size;
                    const int id_shift_left = i_loc_shift * kblock_size;
                    const int i_ub = min(D1, iblock_size - i_loc_shift);

                    for (int jd = 0; jd < jd_ub; ++jd)
                    {
                        const int j_loc_shift = jd * D3;
                        const int j_ub = min(D3, jblock_size - j_loc_shift);
                        const int gen_diff_val = i_ub * j_ub;
                        double* gen_drow = gen_block +
                            id_shift_gen + j_loc_shift * i_ub;

                        for (int kd = 0; kd < kd_ub; ++kd)
                        {
                            const int k_loc_shift = kd * D2;
                            const int k_ub = min(D2, kblock_size - k_loc_shift);
                            const double* left_drow  = left_block +
                                id_shift_left + k_loc_shift * i_ub;
                            const double* right_dblock = right_block +
                                k_loc_shift * jblock_size + j_loc_shift * k_ub;

                            // element level
                            for (int i = 0; i < i_ub; ++i)
                            {
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    double sum = 0.0;
                                    for (int k = 0; k < k_ub; ++k)
                                    {
                                        sum += left_drow[k] * right_dblock[k*j_ub + j];
                                    }
                                    gen_drow[j] += sum;
                                }

                                gen_drow += j_ub;
                                left_drow += k_ub;
                            }

                            gen_drow -= gen_diff_val;
                        }
                    }
                }
            }
        }
    }
}

void matrix_multiplication_block(double* gen_matrix,
    const double* left_matrix, const double* right_matrix,
    const int N1, const int N2, const int N3,
    const int B1, const int B2, const int B3)
{
    const int NB1 = (int)ceil(1.0 * N1 / B1);
    const int NB2 = (int)ceil(1.0 * N2 / B2);
    const int NB3 = (int)ceil(1.0 * N3 / B3);

    // block level
#pragma omp parallel for schedule(dynamic, 1)
    for (int ib = 0; ib < NB1; ++ib)
    {
        const int i_shift_gen = ib*B1*N3;
        const int i_shift_left = ib*B1*N2;
        const int i_ub = min(B1, N1 - ib*B1);

        for (int jb = 0; jb < NB3; ++jb)
        {
            double* gen_block_begin = gen_matrix + i_shift_gen + jb*i_ub*B3;
            const int j_ub = min(B3, N3 - jb*B3);

            for (int kb = 0; kb < NB2; ++kb)
            {
                const int k_ub = min(B2, N2 - kb*B2);
                const double* left_block_begin = left_matrix + i_shift_left + kb*i_ub*B2;
                const double* right_block_begin = right_matrix + kb*B2*N3 + jb*k_ub*B3;

                // element level
                for (int i = 0; i < i_ub; ++i)
                {
                    const double* left_row = left_block_begin + i * k_ub;
                    double* gen_row = gen_block_begin + i * j_ub;
                    for (int j = 0; j < j_ub; ++j)
                    {
                        double sum = 0.0;
                        for (int k = 0; k < k_ub; ++k)
                        {
                            sum += left_row[k] * right_block_begin[k*j_ub + j];
                        }
                        gen_row[j] += sum;
                    }
                }
            }
        }
    }
}

void matrix_multiplication_double_tiled(double* gen_matrix,
    const double* left_matrix, const double* right_matrix,
    const int N1, const int N2, const int N3,
    const int B1, const int B2, const int B3,
    const int D1, const int D2, const int D3)
{
    const int NB1 = (int)ceil(1.0 * N1 / B1);
    const int NB2 = (int)ceil(1.0 * N2 / B2);
    const int NB3 = (int)ceil(1.0 * N3 / B3);
    const int ND1 = (int)ceil(1.0 * B1 / D1);
    const int ND2 = (int)ceil(1.0 * B2 / D2);
    const int ND3 = (int)ceil(1.0 * B3 / D3);

    // big blocks level
#pragma omp parallel for schedule(dynamic, 1)
    for (int ib = 0; ib < NB1; ++ib)
    {
        const int ib_min = min((ib + 1)*B1, N1);
        for (int jb = 0; jb < NB3; ++jb)
        {
            const int jb_min = min((jb + 1)*B3, N3);
            for (int kb = 0; kb < NB2; ++kb)
            {
                const int kb_min = min((kb + 1)*B2, N2);
                // small blocks level
                for (int id = 0; id < ND1; ++id)
                {
                    const int lb1 = ib*B1 + id*D1;
                    const int ub1 = min(ib*B1 + (id + 1)*D1, ib_min);

                    for (int jd = 0; jd < ND3; ++jd)
                    {
                        const int lb3 = jb*B3 + jd*D3;
                        const int ub3 = min(jb*B3 + (jd + 1)*D3, jb_min);

                        for (int kd = 0; kd < ND2; ++kd)
                        {
                            const int lb2 = kb*B2 + kd*D2;
                            const int ub2 = min(kb*B2 + (kd + 1)*D2, kb_min);

                            // elements level
                            for (int i = lb1; i < ub1; ++i)
                            {
                                const double* left_line = left_matrix + i * N2;
                                double* gen_line = gen_matrix + i * N3;
                                for (int j = lb3; j < ub3; ++j)
                                {
                                    double sum = 0.0;
                                    for (int k = lb2; k < ub2; ++k)
                                    {
                                        sum += left_line[k] * right_matrix[k * N3 + j];
                                    }
                                    gen_line[j] += sum;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void matrix_multiplication_tiled(double* gen_matrix,
    const double* left_matrix, const double* right_matrix,
    const int N1, const int N2, const int N3,
    const int B1, const int B2, const int B3)
{
    const int NB1 = (int)ceil(1.0 * N1 / B1);
    const int NB2 = (int)ceil(1.0 * N2 / B2);
    const int NB3 = (int)ceil(1.0 * N3 / B3);

    // block level
#pragma omp parallel for schedule(dynamic, 1)
    for (int ib = 0; ib < NB1; ++ib)
    {
        const int i_lb = ib * B1;
        const int i_ub = min((ib + 1) * B1, N1);
        for (int jb = 0; jb < NB3; ++jb)
        {
            const int j_lb = jb * B3;
            const int j_ub = min((jb + 1) * B3, N3);
            for (int kb = 0; kb < NB2; ++kb)
            {
                const int k_lb = kb * B2;
                const int k_ub = min((kb + 1) * B2, N2);
                // element level
                for (int i = i_lb; i < i_ub; ++i)
                {
                    const double* left_line = left_matrix + i * N2;
                    double* gen_line = gen_matrix + i * N3;
                    for (int j = j_lb; j < j_ub; ++j)
                    {
                        double sum = 0.0;
                        for (int k = k_lb; k < k_ub; ++k)
                        {
                            sum += left_line[k] * right_matrix[k * N3 + j];
                        }
                        gen_line[j] += sum;
                    }
                }
            }
        }
    }
}


void matrix_multiplication_standard(double* gen_matrix,
    const double* src_left_matrix,
    const double* src_right_matrix,
    const int N1, const int N2, const int N3)
{
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < N1; ++i)
    {
        const double* left_line = src_left_matrix + i * N2;
        double* gen_line = gen_matrix + i * N3;
        for (int j = 0; j < N3; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < N2; ++k)
            {
                sum += left_line[k] * src_right_matrix[k * N3 + j];
            }
            gen_line[j] = sum;
        }
    }
}

#ifdef __cplusplus
    }
#endif
