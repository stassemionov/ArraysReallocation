#include "multiplication.h"

#include <algorithm>

using std::min;

void block_matrix_multiplication_double_tiled(double* gen_matrix,
    const double* left_matrix,
    const double* right_matrix,
    const TaskClass& left_mat_data,
    const TaskClass& right_mat_data)
{
    const int N1 = left_mat_data.getDataRef().M_ROWS;
    const int N2 = left_mat_data.getDataRef().M_COLS;
    const int N3 = right_mat_data.getDataRef().M_COLS;

    const int NB1 = left_mat_data.getDataRef().M_BLOCK_ROWS;
    const int NB2 = left_mat_data.getDataRef().M_BLOCK_COLS;
    const int NB3 = right_mat_data.getDataRef().M_BLOCK_COLS;

    const int B1 = left_mat_data.getDataRef().B_ROWS;
    const int B2 = left_mat_data.getDataRef().B_COLS;
    const int B3 = right_mat_data.getDataRef().B_COLS;

    const int D1 = left_mat_data.getDataRef().D_ROWS;
    const int D2 = left_mat_data.getDataRef().D_COLS;
    const int D3 = right_mat_data.getDataRef().D_COLS;

    const int ND1 = static_cast<int>(ceil(1.0 * B1 / D1));
    const int ND2 = static_cast<int>(ceil(1.0 * B2 / D2));
    const int ND3 = static_cast<int>(ceil(1.0 * B3 / D3));

    TaskClass gen_mat_data;
    gen_mat_data.makeData(N1, N3, B1, B3, D1, D3);

    double sum;
    // big block level
    for (int ib = 0; ib < NB1; ++ib)
    {
        const int i_shift_gen = ib*B1*N3;
        const int i_shift_left = ib*B1*N2;

        const int rb1 = min(B1, N1 - ib*B1);
        for (int jb = 0; jb < NB3; ++jb)
        {
            const int rb3 = min(B3, N3 - jb*B3);
            
            double* gen_block_begin = gen_matrix + i_shift_gen + jb*rb1*B3; // begin of big block
            for (int kb = 0; kb < NB2; ++kb)
            {
                const int rb2 = min(B2, N2 - kb*B2);
                
                const double* left_block_begin = left_matrix + i_shift_left + kb*rb1*B2;
                const double* right_block_begin = right_matrix + kb*B2*N3 + jb*rb2*B3;

                const int rnd1 = min(ND1, static_cast<int>(ceil(1.0 * rb1 / D1)));
                const int rnd2 = min(ND2, static_cast<int>(ceil(1.0 * rb2 / D2)));
                const int rnd3 = min(ND3, static_cast<int>(ceil(1.0 * rb3 / D3)));
                // small block level
                for (int id = 0; id < rnd1; ++id)
                {
                    const int id_shift_gen = id*D1*rb3;
                    const int id_shift_left = id*D1*rb2;

                    const int rd1 = min(D1, rb1 - id*D1);
                    for (int jd = 0; jd < rnd3; ++jd)
                    {
                        const int rd3 = min(D3, rb3 - jd*D3);

                        double* gen_double_block_begin = gen_block_begin + id_shift_gen + jd*rd1*D3;
                        for (int kd = 0; kd < rnd2; ++kd)
                        {
                            const int rd2 = min(D2, rb2 - kd*D2);
                            
                            const double* left_double_block_begin  = left_block_begin + id_shift_left + kd*rd1*D2;
                            const double* right_double_block_begin = right_block_begin + kd*D2*rb3 + jd*rd2*D3;
                            // element level
                            for (int i = 0; i < rd1; ++i)
                            {
                                const double* left_line = left_double_block_begin + i*rd2;
                                double* gen_line = gen_double_block_begin + i*rd3;
                                for (int j = 0; j < rd3; ++j)
                                {
                                    sum = gen_line[j];
                                    for (int k = 0; k < rd2; ++k)
                                    {
                                        sum += left_line[k] * right_double_block_begin[k*rd3 + j];
                                    }
                                    gen_line[j] = sum;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void block_matrix_multiplication_tiled(double* gen_matrix,
    const double* left_matrix,
    const double* right_matrix,
    const TaskClass& left_mat_data,
    const TaskClass& right_mat_data)
{
    const int N1 = left_mat_data.getDataRef().M_ROWS;
    const int N2 = left_mat_data.getDataRef().M_COLS;
    const int N3 = right_mat_data.getDataRef().M_COLS;

    const int NB1 = left_mat_data.getDataRef().M_BLOCK_ROWS;
    const int NB2 = left_mat_data.getDataRef().M_BLOCK_COLS;
    const int NB3 = right_mat_data.getDataRef().M_BLOCK_COLS;

    const int B1 = left_mat_data.getDataRef().B_ROWS;
    const int B2 = left_mat_data.getDataRef().B_COLS;
    const int B3 = right_mat_data.getDataRef().B_COLS;

    TaskClass gen_mat_data;
    gen_mat_data.makeData(N1, N3, B1, B3);

    double sum;
    // block level
    for (int ib = 0; ib < NB1; ++ib)
    {
        const int i_shift_gen = ib*B1*N3;
        const int i_shift_left = ib*B1*N2;

        const int rb1 = min(B1, N1 - ib*B1);
        for (int jb = 0; jb < NB3; ++jb)
        {
            const int rb3 = min(B3, N3 - jb*B3);
            double* gen_block_begin = gen_matrix + i_shift_gen + jb*rb1*B3;
            for (int kb = 0; kb < NB2; ++kb)
            {
                const int rb2 = min(B2, N2 - kb*B2);
                // element level
                const double* left_block_begin = left_matrix + i_shift_left + kb*rb1*B2;
                const double* right_block_begin = right_matrix + kb*B2*N3 + jb*rb2*B3;
                
                for (int i = 0; i < rb1; ++i)
                {
                    const double* left_line = left_block_begin + i*rb2;
                    double* gen_line = gen_block_begin + i*rb3;
                    for (int j = 0; j < rb3; ++j)
                    {
                        sum = gen_line[j];
                        for (int k = 0; k < rb2; ++k)
                        {
                            sum += left_line[k] * right_block_begin[k*rb3 + j];
                        }
                        gen_line[j] = sum;
                    }
                }
            }
        }
    }
}

void matrix_multiplication_double_tiled(double*          gen_matrix,
                                        const double*    left_matrix,
                                        const double*    right_matrix,
                                        const TaskClass& left_mat_data,
                                        const TaskClass& right_mat_data)
{
    const int N1 = left_mat_data.getDataRef().M_ROWS;
    const int N2 = left_mat_data.getDataRef().M_COLS;
    const int N3 = right_mat_data.getDataRef().M_COLS;

    const int NB1 = left_mat_data.getDataRef().M_BLOCK_ROWS;
    const int NB2 = left_mat_data.getDataRef().M_BLOCK_COLS;
    const int NB3 = right_mat_data.getDataRef().M_BLOCK_COLS;

    const int B1 = left_mat_data.getDataRef().B_ROWS;
    const int B2 = left_mat_data.getDataRef().B_COLS;
    const int B3 = right_mat_data.getDataRef().B_COLS;

    const int D1 = left_mat_data.getDataRef().D_ROWS;
    const int D2 = left_mat_data.getDataRef().D_COLS;
    const int D3 = right_mat_data.getDataRef().D_COLS;

    const int ND1 = static_cast<int>(ceil(1.0 * B1 / D1));
    const int ND2 = static_cast<int>(ceil(1.0 * B2 / D2));
    const int ND3 = static_cast<int>(ceil(1.0 * B3 / D3));

    double sum;
    // big blocks level
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
                                    sum = gen_line[j];
                                    for (int k = lb2; k < ub2; ++k)
                                    {
                                        sum += left_line[k] * right_matrix[k * N3 + j];
                                    }
                                    gen_line[j] = sum;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void matrix_multiplication_tiled(double*          gen_matrix,
                                 const double*    left_matrix,
                                 const double*    right_matrix,
                                 const TaskClass& left_mat_data,
                                 const TaskClass& right_mat_data)
{
    const int& N1 = left_mat_data.getDataRef().M_ROWS;
    const int& N2 = left_mat_data.getDataRef().M_COLS;
    const int& N3 = right_mat_data.getDataRef().M_COLS;

    const int& NB1 = left_mat_data.getDataRef().M_BLOCK_ROWS;
    const int& NB2 = left_mat_data.getDataRef().M_BLOCK_COLS;
    const int& NB3 = right_mat_data.getDataRef().M_BLOCK_COLS;

    const int& B1 = left_mat_data.getDataRef().B_ROWS;
    const int& B2 = left_mat_data.getDataRef().B_COLS;
    const int& B3 = right_mat_data.getDataRef().B_COLS;

    double sum;
    // block level
    for (int ib = 0; ib < NB1; ++ib)
    {
        for (int jb = 0; jb < NB3; ++jb)
        {
            for (int kb = 0; kb < NB2; ++kb)
            {
                // element level
                for (int i = ib*B1; i < min((ib+1)*B1,N1); ++i)
                {
                    const double* left_line = left_matrix + i * N2;
                    double* gen_line = gen_matrix + i * N3;
                    for (int j = jb*B3; j < min((jb + 1)*B3, N3); ++j)
                    {
                        sum = gen_line[j];
                        for (int k = kb*B2; k <  min((kb + 1)*B2, N2); ++k)
                        {
                            sum += left_line[k] * right_matrix[k * N3 + j];
                        }
                        gen_line[j] = sum;
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
    double sum;
    for (int i = 0; i < N1; ++i)
    {
        const double* left_line = src_left_matrix + i * N2;
        double* gen_line = gen_matrix + i * N3;
        for (int j = 0; j < N3; ++j)
        {
            sum = 0;
            for (int k = 0; k < N2; ++k)
            {
                sum += left_line[k] * src_right_matrix[k * N3 + j];
            }
            gen_line[j] = sum;
        }
    }
}
