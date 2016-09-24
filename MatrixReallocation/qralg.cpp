#include "qralg.h"
#include "service.h"

#include <algorithm>

using std::min;
using std::max;
using std::swap;

double* QR_WY_tiled(double* A, const TaskData& parameters)
{
    // Main parameters of task
    const int N = parameters.M_ROWS;
    const int b1 = parameters.B_ROWS;
    const int b2 = parameters.B_COLS;
    
    // Counts of blocks in matrix
    const int block_count_in_col = static_cast<int>(ceil(1.0 * N / b1));
    const int block_count_in_row = static_cast<int>(ceil(1.0 * N / b2));
    
    double* work_memory = reinterpret_cast<double*>(
            calloc(N + b2 + 2*N*b2 + 2*N*N, sizeof(double)));

    double* v =       new(work_memory)   double[N];
    double* w =       new(v + N)         double[b2];
    double* W =       new(w + b2)        double[N*b2];
    double* Y =       new(W + N*b2)      double[N*b2];
    double* WY =      new(Y + N*b2)      double[N*N];
    double* Abuf =    new(WY + N*N)      double[N*N];

    // Indices of current block column and stripe
    int curr_i_block_ind = -1;
    int curr_j_block_ind = -1;

    // 'lambda' is index of first column of current block-column
    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // Indices of big block, that contains element (lambda,lambda).
        // Remark: row-index value can increase with passing over diagonal
        //         when 'b1' != 'b2'.
        curr_i_block_ind = lambda / b1;
        ++curr_j_block_ind;
        // Vertical shift to diagonal block
        const int d1_shift = curr_i_block_ind * b1;
        // Index of first column of next block-column
        const int t = min(lambda + b2, N);
        // Real values of height of current block-stripe/column
        const int d2 = t - lambda;
        // Count of big blocks from the right of current big block in its block-stripe
        const int row_b_count = block_count_in_row - curr_j_block_ind - 1;
        // Count of big blocks below and including current big block in its block-column
        const int col_b_count = block_count_in_col - curr_i_block_ind;

        // Выполнение преобразований над текущим блоком
        for (int it = lambda; it < t; ++it)
        {
            const int loc_it = it - lambda;
            const int w_diff_val = t - it;
            const int A_diff_val = N - w_diff_val;
            // Индекс большого блока, содержащего 'lambda+it'-ый диагональный элемент,
            // считая от блока, содержащего 'lambda'-ый элемент
            const int b_begin = (it - d1_shift) / b1;

            memset(w, 0, d2*sizeof(double));

            double norm = 0.0;
            double scalar = 1.0;
            // * Вычисление вектора Хаусхолдера
            const double* A_shifted_HV = A + it*N + it;
            double* v_shifted_HV = v + it;
            for (int j = 0; j < N - it; ++j)
            {
                const double buf(*A_shifted_HV);
                *v_shifted_HV++ = buf;
                norm += buf * buf;
                A_shifted_HV += N;
            }
            const double A_diag_el = A[it*N + it];
            const double diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + diag_el_sign * sqrt(norm);
            // используется нормировка, т.ч. v[it] = 1
            for (int j = it + 1; j < N; ++j)
            {
                double buf(v[j] /= beta);
                scalar += buf * buf;
            }
            v[it] = 1.0;
            beta = -2.0 / scalar;

            // * Применение преобразования

            // Вычисление вектора w = beta * (A(:,lambda:t-1)^t) * v .
            // начинаем с it, т.к. it-ая матрица Хаусхолдера действует на матрицу A(it:,it:)
            const int j_iter_count_TRANSF = t - it;
            for (int ib = b_begin; ib < col_b_count; ++ib)
            {
                const int ib_shift = d1_shift + ib*b1;
                const int lb = max(it, ib_shift);
                const int ub = min(ib_shift + b1, N);
                const double* A_row_shifted = A + lb*N + it;
                const double* v_shifted = v + lb;
                double* w_shifted = w + loc_it;

                for (int i = 0; i < ub - lb; ++i)
                {
                    const double vi(*v_shifted++);
                    for (int j = 0; j < j_iter_count_TRANSF; ++j)
                    {
                        (*w_shifted++) += (*A_row_shifted++) * vi;
                    }
                    A_row_shifted += A_diff_val;
                    w_shifted -= w_diff_val;
                }
            }
            for (int k = loc_it; k < d2; ++k)
            {
                w[k] *= beta;
            }

            // Transformation of current block-column:
            // A(it:N, loc_it:d2-1) += v[it:N] * (w[loc_it:d2-1]^t)
            for (int ib = b_begin; ib < col_b_count; ++ib)
            {
                const int mb_shift = d1_shift + ib * b1;
                const int i_lb = max(it, mb_shift);
                const int i_ub = min(mb_shift + b1, N);
                const double* v_shifted = v + i_lb;
                const double* w_shifted = w + loc_it;
                double* A_row_shifted = A + i_lb*N + it;

                for (int i = 0; i < i_ub - i_lb; ++i)
                {
                    const double vi(*v_shifted++);
                    for (int j = 0; j < j_iter_count_TRANSF; ++j)
                    {
                        (*A_row_shifted++) += vi * (*w_shifted++);
                    }
                    A_row_shifted += A_diff_val;
                    w_shifted -= w_diff_val;
                }
            }

            // * Заполнение поддиагнальной части it-го столбца матрицы A
            //   информативной частью вектора Хаусхолдера
            double* A_shifted_fill = A + (it + 1)*N + it;
            const double* v_shifted_fill = v + it + 1;
            for (int j = it + 1; j < N; ++j)
            {
                (*A_shifted_fill) = (*v_shifted_fill++);
                A_shifted_fill += N;
            }
        }// FOR

         // * Вычисление WY-представления произведения матриц Хаусхолдера
        const int d = t - lambda;
        for (int it = 0; it < d; ++it)
        {
            const int shift = lambda + it;
            const int b_begin = (shift - d1_shift) / b1;

            memset(v + d1_shift, 0, (N - d1_shift)*sizeof(double));

            double scalar = 1, buf;
            // * Вычисление нормы и ск. произведения вектора Хаусхолдера
            const double* A_shifted_HV = A + (shift + 1)*N + shift;
            double* v_shifted_HV = v + shift + 1;
            for (int j = shift + 1; j < N; ++j)
            {
                buf = *A_shifted_HV;
                scalar += buf * buf;
                (*v_shifted_HV++) = buf;
                A_shifted_HV += N;
            }
            v[shift] = 1.0;
            double beta = -2.0 / scalar;

            // it - число накопленных столбцов в матрицах W и Y (= номеру дозаписываемого столбца)
            if (it == 0)
            {
                // Начальное заполнение первого столбца
                // замечание: первые lambda строк матриц W и Y - нулевые,
                //            матрицы имеют ступенчатый вид.
                double* Y_shifted_init = Y + lambda*b2;
                double* W_shifted_init = W + lambda*b2;
                const double* v_shifted_init = v + lambda;
                for (int i = 0; i < N - lambda; ++i)
                {
                    const double buf_v(*v_shifted_init++);
                    *Y_shifted_init = buf_v;
                    *W_shifted_init = beta * buf_v;
                    Y_shifted_init += b2;
                    W_shifted_init += b2;
                }
            }
            else
            {
                memset(w, 0, d * sizeof(double));

                const int W_Y_diff_val = b2 - it;

                // Вычисление произведения (Y^t) * v
                double* wp_Yv = w;
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    const int i_lb = max(shift, d1_shift + ib * b1);
                    const int i_ub = min(d1_shift + (ib + 1) * b1, N);
                    const double* Y_row = Y + i_lb * b2;
                    const double* v_shifted = v + i_lb;

                    for (int i = 0; i < i_ub - i_lb; ++i)
                    {
                        const double vi(*v_shifted++);
                        for (int j = 0; j < it; ++j)
                        {
                            (*wp_Yv++) += (*Y_row++) * vi;
                        }
                        wp_Yv -= it;
                        Y_row += W_Y_diff_val;
                    }
                }

                // Вычисление произведения W * ((Y^t) * v) <=> W * w
                const double* wp_Ww = w;
                for (int ib = 0; ib < col_b_count; ++ib)
                {
                    const int i_lb = d1_shift + ib * b1;
                    const int i_ub = min(i_lb + b1, N);
                    double* W_row = W + i_lb * b2;
                    double* Y_row_shifted = Y + i_lb * b2 + it;
                    double* v_shifted = v + i_lb;

                    for (int i = 0; i < i_ub - i_lb; ++i)
                    {
                        const double vi = *v_shifted++;
                        double sum(0.0);
                        for (int j = 0; j < it; ++j)
                        {
                            sum += (*W_row++) * (*wp_Ww++);
                        }
                        *W_row = beta * (vi + sum);
                        *Y_row_shifted = vi;

                        wp_Ww -= it;
                        W_row += W_Y_diff_val;
                        Y_row_shifted += b2;
                    }
                }
            }
        }

        // * Преобразование остальной части матрицы A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        // Multiplication Y * W^t
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib * b1;
            const int i_lb = max(lambda, i_mb_shift);
            const int i_ub = min(i_mb_shift + b1, N);
            const int i_iter_count = i_ub - i_lb;
            const int Y_diff_val = i_iter_count * b2;
            const int i_lb_mul_N = i_lb * N;
            const double* Y_row = Y + i_lb * b2;

            for (int jb = 0; jb < col_b_count; ++jb)
            {
                const int j_mb_shift = d1_shift + jb * b1;
                const int j_lb = max(lambda, j_mb_shift);
                const int j_ub = min(j_mb_shift + b1, N);
                const int j_iter_count = j_ub - j_lb;
                const int WY_diff_val = N - j_iter_count;
                const int W_diff_val = j_iter_count * b2;
                double* WY_row_shifted = WY + i_lb_mul_N + j_lb;
                const double* W_row = W + j_lb * b2;

                for (int i = 0; i < i_iter_count; ++i)
                {
                    const int kbound = min(d, i_lb + i - lambda + 1);
                    const int Wrow_diff_val = b2 - kbound;
                    for (int j = 0; j < j_iter_count; ++j)
                    {
                        double sum(0.0);
                        for (int k = 0; k < kbound; ++k)
                        {
                            sum += (*Y_row++) * (*W_row++);
                        }
                        (*WY_row_shifted++) = sum;

                        Y_row -= kbound;
                        W_row += Wrow_diff_val;
                    }
                    Y_row += b2;
                    W_row -= W_diff_val;
                    WY_row_shifted += WY_diff_val;
                }
                Y_row -= Y_diff_val;
            }
        }

        const double* A_shifted_pre_copy = A + lambda * N + t;
        double* Abuf_shifted_pre_copy = Abuf + lambda * N + t;
        const int copy_length = (N - t) * sizeof(double);
        for (int i = lambda; i < N; ++i)
        {
            memcpy(Abuf_shifted_pre_copy, A_shifted_pre_copy, copy_length);
            Abuf_shifted_pre_copy += N;
            A_shifted_pre_copy += N;
        }

        // Умножение A(lamda:N-1, t:N-1) += (Y * W^t) * A(lamda:N-1, t:N-1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_lb = lambda + ib * b1;
            const int i_ub = min(i_lb + b1, N);
            const int i_iter_count = i_ub - i_lb;
            const int Abuf_diff_val = i_iter_count * N;
            
            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_lb = t + jb * b2;
                const int j_ub = min(j_lb + b2, N);
                const int j_iter_count = j_ub - j_lb;
                const int Abuf_row_diff_val = N - j_iter_count;
                double* Abuf_row_shifted = Abuf + i_lb * N + j_lb;

                for (int kb = 0; kb < col_b_count; ++kb)
                {
                    const int k_lb = lambda + kb * b1;
                    const int k_ub = min(k_lb + b1, N);
                    const int k_iter_count = k_ub - k_lb;
                    const int A_row_diff_val = k_iter_count * N - 1;
                    const double* WY_row_shifted = WY + i_lb * N + k_lb;
                    const double* A_row_shifted = A + k_lb * N + j_lb;

                    for (int i = 0; i < i_iter_count; ++i)
                    {
                        for (int j = 0; j < j_iter_count; ++j)
                        {
                            double sum(*Abuf_row_shifted);
                            for (int k = 0; k < k_iter_count; ++k)
                            {
                                sum += (*WY_row_shifted++) * (*A_row_shifted);
                                A_row_shifted += N;
                            }
                            (*Abuf_row_shifted++) = sum;

                            WY_row_shifted -= k_iter_count;
                            A_row_shifted -= A_row_diff_val;
                        }
                        Abuf_row_shifted += Abuf_row_diff_val;
                        WY_row_shifted += N;
                        A_row_shifted -= j_iter_count;
                    }
                    Abuf_row_shifted -= Abuf_diff_val;
                }
            }
        }
                
        // копирование преобразованной части матрицы A
        double* A_shifted_post_copy = A + lambda * N + t;
        const double* Abuf_shifted_post_copy = Abuf + lambda * N + t;
        for (int i = lambda; i < N; ++i)
        {
            memcpy(A_shifted_post_copy, Abuf_shifted_post_copy, copy_length);
            Abuf_shifted_post_copy += N;
            A_shifted_post_copy += N;
        }

    }// WHILE

    free(work_memory);

    return A;
}

double* QR_WY_double_tiled(double* A, const TaskData& parameters)
{
    // Main parameters of data layout
    const int N = parameters.M_ROWS;
    const int b1 = parameters.B_ROWS;
    const int b2 = parameters.B_COLS;
    const int db1 = parameters.D_ROWS;
    const int db2 = parameters.D_COLS;

    const int bsizes_ratio = static_cast<int>(ceil(1.0 * b2 / b1));
    // Big blocks counts in matrix
    const int block_count_in_col = static_cast<int>(ceil(1.0 * N / b1));
    const int block_count_in_row = static_cast<int>(ceil(1.0 * N / b2));
    // Small blocks counts in big blocks
    const int dblock_count_in_bcol = static_cast<int>(ceil(1.0 * b1 / db1));
    const int dblock_count_in_brow = static_cast<int>(ceil(1.0 * b2 / db2));
    const int dblock_count_in_diff_bcol = (N % b1 == 0) ?
        dblock_count_in_bcol :
        static_cast<int>(ceil(1.0 * (N % b1) / db1));
    const int dblock_count_in_diff_brow = (N % b2 == 0) ?
        dblock_count_in_brow :
        static_cast<int>(ceil(1.0 * (N % b2) / db2));

    double* work_memory = reinterpret_cast<double*>(
            calloc(N + b2 + 2*N*b2 + 2*N*N + db1, sizeof(double)));

    double* v =       new(work_memory)   double[N];
    double* w =       new(v + N)         double[b2];
    double* sum_vec = new(w + b2)        double[db1];
    double* W =       new(sum_vec + db1) double[N*b2];
    double* Y =       new(W + N*b2)      double[N*b2];
    double* WY =      new(Y + N*b2)      double[N*N];
    double* Abuf =    new(WY + N*N)      double[N*N];

    int curr_i_block_ind = -1;
    int curr_j_block_ind = -1;

    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // Indices of big block, that contains element (lambda,lambda).
        // Remark: row-index value can increase with passing over diagonal
        //         when 'b1' != 'b2'.
        curr_i_block_ind = lambda / b1;
        ++curr_j_block_ind;

        // Index of first column of next block-column
        const int t = min(lambda + b2, N);
        // Count of big blocks from the right of current big block in its block-stripe
        const int row_b_count = block_count_in_row - curr_j_block_ind - 1;
        // Count of big blocks below and including current big block in its block-column
        const int col_b_count = block_count_in_col - curr_i_block_ind;
        // Horizontal and vertical shifts to diagonal block
        const int d1_shift = curr_i_block_ind * b1;
        const int d2_shift = curr_j_block_ind * b2;
        // Real values of height/width of current block-stripe/column
        const int d2 = t - lambda;
        // Count of small blocks in current big block in row-direction
        const int dblock_count_in_curr_block_row =
            (curr_j_block_ind == block_count_in_row - 1) ?
            dblock_count_in_diff_brow : dblock_count_in_brow;

        // * Execution of current block-column transformations...
        for (int it = 0; it < d2; ++it)
        {
            memset(w, 0, d2 * sizeof(double));
            memset(v + lambda, 0, it * sizeof(double));

            // Shift for 'lambda+it'-th element of matrix column
            const int full_shift = lambda + it;
            // Индекс большого блока, содержащего 'lambda+it'-ый диагональный элемент,
            // считая от блока, содержащего 'lambda'-ый элемент
            const int b_begin = (full_shift - d1_shift) / b1;
            // Index of small-block-column,
            // which contains 'lambda+it'-th column
            const int dblock_col_index = it / db2;

            // * Householder vector computation...
            double norm = 0.0;
            double scalar = 1.0;
            const double* A_shifted_HV = A + full_shift * N + full_shift;
            double* v_shifted_HV = v + full_shift;
            for (int i = 0; i < N - full_shift; ++i)
            {
                double buf(*A_shifted_HV);
                v_shifted_HV[i] = buf;
                norm += buf * buf;
                A_shifted_HV += N;
            }
            const double A_diag_el = A[full_shift*N + full_shift];
            const double A_diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + A_diag_el_sign * sqrt(norm);
            // используется нормировка, т.ч. v[it] = 1
            double* v_shifted_NORM = v + full_shift + 1;
            for (int j = 0; j < N-full_shift-1; ++j)
            {
                double buf(v_shifted_NORM[j] /= beta);
                scalar += buf * buf;
            }
            v[full_shift] = 1.0;
            beta = -2.0 / scalar;

            // * Transformation applying...

            // Computation of vector w = beta * (A(:,lambda:t-1)^t) * v .
            // Starting with 'full_shift'-th row and column of matrix,
            // because 'full_shift'-th Householder matrix really operates
            // with A(full_shift:*, full_shift:t-1).
            // Pattern of data access is as follows:
            // 1. Fix big block (ib)
            // 2. Fix small block stripe in fixed big block (id)
            // 3. Fix small block in fixed stripe (jd)
            // 4. Fix row in fixed small block (i)
            // 5. Produce passing of fixed row (j)
            for (int ib = b_begin; ib < col_b_count; ++ib)
            {
                const int mb_shift = d1_shift + ib * b1;
                const int block_lb = max(full_shift, mb_shift) - mb_shift;
                const int block_ub = min(mb_shift + b1, N);
                const int id_lb = (block_lb == 0) ? 0 : block_lb / db1;
                const int id_ub = (ib == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;

                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int i_bd_shift = mb_shift + id * db1;
                    const int i_lb = max(full_shift, i_bd_shift) - i_bd_shift;
                    const int i_ub = min(db1, block_ub - i_bd_shift);
                    const int i_iter_count = i_ub - i_lb;
                    const double* A_dstripe = A + i_bd_shift * N;
                    const double* v_id_shifted = v + i_bd_shift + i_lb;

                    for (int jd = dblock_col_index; jd < dblock_count_in_curr_block_row; ++jd)
                    {
                        const int loc_drow_shift = jd * db2;
                        const int j_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int j_ub = min(d2 - loc_drow_shift, db2);
                        const int j_iter_count = j_ub - j_lb;
                        const int A_dif_val = N - j_iter_count;
                        double* w_shifted = w + loc_drow_shift + j_lb;
                        const double* A_row = (A_dstripe + i_lb*N) +
                            lambda + loc_drow_shift + j_lb;

                        for (int i = 0; i < i_iter_count; ++i)
                        {
                            const double vi = *v_id_shifted++;
                            for (int j = j_lb; j < j_ub; ++j)
                            {
                                (*w_shifted++) += (*A_row++) * vi;
                            }

                            w_shifted -= j_iter_count;
                            A_row += A_dif_val;
                        }
                        v_id_shifted -= i_iter_count;
                    }
                }
            }
            for (int k = it; k < d2; ++k)
            {
                w[k] *= beta;
            }

            // Transformation of current block-column:
            // A(fullshift:N, it:t-1) += v[fullshift:N] * (w[it:t-1]^t)
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                const int mb_shift = d1_shift + jb * b1;
                const int block_lb = max(full_shift, mb_shift) - mb_shift;
                const int block_ub = min(mb_shift + b1, N);
                const int jd_lb = (block_lb == 0) ? 0 : block_lb / db1;
                const int jd_ub = (jb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;

                for (int jd = jd_lb; jd < jd_ub; ++jd)
                {
                    const int bd_shift = mb_shift + jd * db1;
                    const int j_lb = max(full_shift, bd_shift) - bd_shift;
                    const int j_ub = min(db1, block_ub - bd_shift);
                    const int j_iter_count = j_ub - j_lb;
                    const int v_diff_val = j_ub - j_lb;
                    const int j_lb_mult_N = j_lb * N;
                    double* A_dstripe = A + bd_shift * N;
                    const double* v_jd_shifted = v + bd_shift + j_lb;

                    for (int kd = dblock_col_index; kd < dblock_count_in_curr_block_row; ++kd)
                    {
                        const int loc_drow_shift = kd * db2;
                        const int k_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int dblock_w = min(db2, d2 - loc_drow_shift);
                        const int k_iter_count = dblock_w - k_lb;
                        const int A_diff_val = N - k_iter_count;
                        const double* w_shifted = w + loc_drow_shift + k_lb;
                        double* A_line = (A_dstripe + j_lb_mult_N) +
                            lambda + loc_drow_shift + k_lb;

                        for (int j = 0; j < j_iter_count; ++j)
                        {
                            const double vj = (*v_jd_shifted++);
                            for (int k = 0; k < k_iter_count; ++k)
                            {
                                (*A_line++) += vj * (*w_shifted++);
                            }
                            A_line += A_diff_val;
                            w_shifted -= k_iter_count;
                        }

                        v_jd_shifted -= v_diff_val;
                    }
                }
            }

            // * Заполнение поддиагнальной части it-го столбца матрицы A
            //   информативной частью вектора Хаусхолдера
            double* A_shifted_fill = A + (full_shift + 1)*N + full_shift;
            const double* v_shifted_fill = v + full_shift + 1;
            for (int j = full_shift + 1; j < N; ++j)
            {
                (*A_shifted_fill) = (*v_shifted_fill++);
                A_shifted_fill += N;
            }
        }// FOR

        // * Computation of WY-representation of Householder matrices production...

        const int set_len = min(bsizes_ratio*b1, N - d1_shift)*b2*sizeof(double);
        memset(W + d1_shift*b2, 0, set_len);
        memset(Y + d1_shift*b2, 0, set_len);
        for (int it = 0; it < d2; ++it)
        {
            memset(w, 0, d2*sizeof(double));
            memset(v + lambda, 0, it*sizeof(double));

            const int full_shift = lambda + it;
            const int b_begin = (full_shift - d1_shift) / b1;
            const int dblock_col_index = it / db2;

            // Вычисление нормы и ск. произведения вектора Хаусхолдера
            double scalar(1.0);
            double* A_shifted_HV = A + (full_shift + 1)*N + full_shift;
            double* v_shifted_HV = v + full_shift + 1;
            for (int j = full_shift + 1; j < N; ++j)
            {
                double buf(*A_shifted_HV);
                scalar += buf * buf;
                (*v_shifted_HV++) = buf;
                A_shifted_HV += N;
            }
            v[full_shift] = 1.0;
            double beta(-2.0 / scalar);

            // 'it' is count of accumulated columns in W and Y matrices
            // (equals index of columns, which is added at current iteration)
            if (it == 0)
            {
                // Initial first column filling.
                // Remark: first 'lambda' of rows in matrices W and Y are filled with zero,
                //           next 'd2' of rows consist lower triangular matrix.
                double* Y_shifted = Y + lambda*b2;
                double* W_shifted = W + lambda*b2;
                const double* v_shifted_WYinit = v + lambda;
                for (int i = lambda; i < N; ++i)
                {
                    double buf_WYinit(*v_shifted_WYinit++);
                    *Y_shifted = buf_WYinit;
                    *W_shifted = beta * buf_WYinit;
                    Y_shifted += b2;
                    W_shifted += b2;
                }
            }
            else
            {
                // Computation of production w = (Y^t) * v .
                // Pattern of data access is as follows:
                // 1. Fix big block (ib)
                // 2. Fix small block stripe in fixed big block (id)
                // 3. Fix small block in fixed stripe (jd)
                // 4. Fix row in fixed small block (i) and element of vector 'v'
                // 5. Produce passing of fixed row of Y (j) and vector 'w',
                //    using fixed element of vector 'v'
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    const int mb_shift = d1_shift + ib * b1;
                    const int block_lb = max(full_shift, mb_shift) - mb_shift;
                    const int block_ub = min(mb_shift + b1, N);
                    const int id_lb = (block_lb == 0) ? 0 : (block_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;

                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int bd_shift = mb_shift + id * db1;
                        const int i_lb = max(full_shift, bd_shift) - bd_shift;
                        const int i_ub = min(db1, block_ub - bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int v_diff_val = i_ub - i_lb;
                        const int i_lb_mult_b2 = i_lb * b2;
                        const double* Ydstripe = Y + bd_shift*b2;
                        const double* v_id_shifted = v + bd_shift + i_lb;

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it - loc_drow_shift, db2);
                            const int Y_diff_val = b2 - j_ub;
                            const double* Y_line = (Ydstripe + i_lb_mult_b2) + loc_drow_shift;
                            double* w_shifted = w + loc_drow_shift;

                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                const double vi = (*v_id_shifted++);
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    (*w_shifted++) += (*Y_line++) * vi;
                                }
                                Y_line += Y_diff_val;
                                w_shifted -= j_ub;
                            }

                            v_id_shifted -= v_diff_val;
                        }
                    }
                }

                // Computation of production W * ((Y^t) * v) <=> W * w .
                // Pattern of data access is as follows:
                // 1. Fix big block (ib)
                // 2. Fix small block stripe in fixed big block (id)
                // 3. Fix small block in fixed stripe (jd)
                // 4. Fix row in fixed small block (i)
                // 5. Produce passing of fixed row of W (j) and vector 'w'
                const int set_len_Ww = db1 * sizeof(double);
                for (int ib = 0; ib < col_b_count; ++ib)
                {
                    const int mb_shift = d1_shift + ib * b1;
                    const int iblock_lb = max(lambda, mb_shift) - mb_shift;
                    const int iblock_ub = min(mb_shift + b1, N);
                    const int id_lb = (iblock_lb == 0) ? 0 : (iblock_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;

                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int bd_shift = mb_shift + id * db1;
                        const int i_lb = max(lambda, bd_shift) - bd_shift;
                        const int i_ub = min(db1, iblock_ub - bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int i_lb_mult_b2 = i_lb * b2;
                        double* Wdstripe = W + bd_shift * b2;
                        double* Ydstripe = Y + (Wdstripe - W);
                        // Vector to contain results of production of W's row and vector w
                        memset(sum_vec, 0, set_len_Ww);

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it, loc_drow_shift + db2) - loc_drow_shift;
                            const int diff_val = b2 - j_ub;
                            const double* w_shifted = w + loc_drow_shift;
                            const double* Wrow = Wdstripe + i_lb_mult_b2 + loc_drow_shift;
                            double* sum_vec_shifted = sum_vec + i_lb;

                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                double sum_loc(0.0);
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    sum_loc += (*Wrow++) * (*w_shifted++);
                                }
                                // Scalar production of W's row and vector w
                                (*sum_vec_shifted++) += sum_loc;

                                Wrow += diff_val;
                                w_shifted -= j_ub;
                            }
                        }
                        // Insertion of computed values into 'it'-th column of Y and W
                        const double* v_shifted_WYadd = v + bd_shift + i_lb;
                        // Pointers to 'it'-th column of Y and W
                        double* Wdblock_add = Wdstripe + i_lb_mult_b2 + it;
                        double* Ydblock_add = Ydstripe + (Wdblock_add - Wdstripe);
                        const double* sum_vec_shifted = sum_vec + i_lb;
                        for (int i = 0; i < i_iter_count; ++i)
                        {
                            double buf(v_shifted_WYadd[i]);
                            *Wdblock_add = beta * (buf + sum_vec_shifted[i]);
                            *Ydblock_add = buf;
                            Wdblock_add += b2;
                            Ydblock_add += b2;
                        }
                    }
                }
            }
        }

        // * Transformation of remaining block-columns of A...
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        double* WY_shifted = WY + d1_shift*N + d1_shift;
        const int set_length = (N - d1_shift)*sizeof(double);
        for (int i = d1_shift; i < N; ++i)
        {
            memset(WY_shifted, 0, set_length);
            WY_shifted += N;
        }

        // Multiplication Y * W^t.
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int iblock_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int iblock_ub = min(i_mb_shift + b1, N);
            const int id_lb = (iblock_lb == 0) ? 0 : (iblock_lb / db1);
            const int id_ub = (ib == col_b_count - 1) ?
                dblock_count_in_diff_bcol :
                dblock_count_in_bcol;

            for (int jb = 0; jb < col_b_count; ++jb)
            {
                const int j_mb_shift = d1_shift + jb*b1;
                const int jblock_lb = max(lambda, j_mb_shift) - j_mb_shift;
                const int jblock_ub = min(j_mb_shift + b1, N);
                const int jd_lb = (jblock_lb == 0) ? 0 : (jblock_lb / db1);
                const int jd_ub = (jb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;

                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int i_bd_shift = i_mb_shift + id * db1;
                    const int i_lb = max(lambda, i_bd_shift) - i_bd_shift;
                    const int i_ub = min(db1, iblock_ub - i_bd_shift);
                    const int i_iter_count = i_ub - i_lb;
                    const int i_lb_mul_b2 = i_lb * b2;
                    const int i_lb_mul_N = i_lb * N;
                    const double* Ydstripe = Y + i_bd_shift*b2;
                    double* WYdstripe = WY + i_bd_shift*N;

                    // Upper bound for 'kd'-level small-block-cycle is defined by
                    // lower triangular form of matrix Y.
                    // If 'id'-th small-block-stripe is crossing matrix diagonal,
                    // then upper bound for small-block-stripe passing equals
                    // index of small block, that is containing diagonal element
                    // with maximum index within 'id'-th small-block-stripe height bounds.
                    const int kd_ub = (i_bd_shift < t) ?
                        ((min(i_bd_shift + db1, iblock_ub) - lambda) / db2 + 1) :
                        dblock_count_in_curr_block_row;

                    for (int jd = jd_lb; jd < jd_ub; ++jd)
                    {
                        const int j_bd_shift = j_mb_shift + jd * db1;
                        const int j_lb = max(lambda, j_bd_shift) - j_bd_shift;
                        const int j_ub = min(db1, jblock_ub - j_bd_shift);
                        const int j_lb_mul_b2 = j_lb * b2;
                        const int WY_diff_val = N - (j_ub - j_lb);
                        const double* Wdstripe = W + j_bd_shift*b2;

                        for (int kd = 0; kd < kd_ub; ++kd)
                        {
                            const int k_loc_shift = kd * db2;
                            const int k_ub = min(d2 - k_loc_shift, db2);
                            const int W_diff_val = b2 - k_ub;

                            // Pointers to first row of small tile of W, Y and WY.
                            // Declared there because maust be updated on each iteration of 'kd'.
                            double* WYrow = (WYdstripe + i_lb_mul_N) + j_bd_shift + j_lb;
                            const double* Yrow = (Ydstripe + i_lb_mul_b2) + k_loc_shift;
                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                const double* Wrow = (Wdstripe + j_lb_mul_b2) + k_loc_shift;
                                for (int j = j_lb; j < j_ub; ++j)
                                {
                                    double sum = 0.0;
                                    for (int k = 0; k < k_ub; ++k)
                                    {
                                        sum += (*Yrow++) * (*Wrow++);
                                    }
                                    (*WYrow++) += sum;

                                    Yrow -= k_ub;
                                    Wrow += W_diff_val;
                                }
                                Yrow += b2;
                                WYrow += WY_diff_val;
                            }
                        }
                    }
                }
            }
        }

        const double* A_shifted = A + d1_shift*N + d2_shift;
        double* Abuf_shifted = Abuf + (A_shifted - A);
        const int copy_length = (N - d2_shift) * sizeof(double);
        for (int i = d1_shift; i < N; ++i)
        {
            memcpy(Abuf_shifted, A_shifted, copy_length);
            A_shifted += N;
            Abuf_shifted += N;
        }

        // Multiplication A(lamda:N-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1) <=>
        //                A(lamda:N-1, t:N-1) += WY * A(lamda:N-1, t:N-1) .
        // Layout pattern: (b1 x b2, d1 x d2) = (b1 x b1, d1 x d1) * (b1 x b2, d1 x d2).
        // Result of multiplication is added to temporal buffer 'Abuf',
        // which is already contain elements of 'A'.
        // It allows to avoid using of elements of matrix 'A',
        // which have results of multiplication, when we need their old values.
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int iblock_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int iblock_ub = min(i_mb_shift + b1, N);
            const int id_lb = (iblock_lb == 0) ? 0 : (iblock_lb / db1);
            const int id_ub = (ib == col_b_count - 1) ?
                dblock_count_in_diff_bcol :
                dblock_count_in_bcol;

            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_mb_shift = t + jb*b2;
                const int jblock_ub = min(j_mb_shift + b2, N);
                const int jd_ub = (jb == row_b_count - 1) ?
                    dblock_count_in_diff_brow :
                    dblock_count_in_brow;

                for (int kb = 0; kb < col_b_count; ++kb)
                {
                    const int k_mb_shift = d1_shift + kb*b1;
                    const int kblock_lb = max(lambda, k_mb_shift) - k_mb_shift;
                    const int kblock_ub = min(k_mb_shift + b1, N);
                    const int kd_lb = (kblock_lb == 0) ? 0 : (kblock_lb / db1);
                    const int kd_ub = (kb == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;

                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int i_bd_shift = i_mb_shift + id * db1;
                        const int i_lb = max(lambda, i_bd_shift) - i_bd_shift;
                        const int i_ub = min(db1, iblock_ub - i_bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int i_lb_mul_N = i_lb * N;
                        double* Abuf_dstripe = Abuf + i_bd_shift*N;
                        const double* WY_dstripe = WY + (Abuf_dstripe - Abuf);

                        for (int jd = 0; jd < jd_ub; ++jd)
                        {
                            const int j_bd_shift = j_mb_shift + jd * db2;
                            const int j_ub = min(j_bd_shift + db2, jblock_ub) - j_bd_shift;
                            const int Abuf_dif_value = N - j_ub;
                            // Points to first row of ['id','jd'] small tile of 'Abuf'
                            double* Abuf_dstripe_shifted = (Abuf_dstripe + i_lb_mul_N) + j_bd_shift;

                            for (int kd = kd_lb; kd < kd_ub; ++kd)
                            {
                                const int k_bd_shift = k_mb_shift + kd * db1;
                                const int k_ub = min(k_bd_shift + db1, kblock_ub) - k_bd_shift;
                                const int A_diff_val = N * k_ub - 1;
                                // Points to 'i_lb'-th row of ['id','kd'] small tile of 'WY'
                                const double* WY_dstripe_shifted = (WY_dstripe + i_lb_mul_N) + k_bd_shift;
                                // Points to first row of ['kd','jd'] small tile of 'A'
                                const double* A_row = (A + k_bd_shift*N) + j_bd_shift;

                                double* Abuf_row = Abuf_dstripe_shifted;
                                for (int i = 0; i < i_iter_count; ++i)
                                {
                                    const double* WY_row = WY_dstripe_shifted;
                                    for (int j = 0; j < j_ub; ++j)
                                    {
                                        double sum(0.0);
                                        for (int k = 0; k < k_ub; ++k)
                                        {
                                            sum += (*WY_row++) * (*A_row);
                                            A_row += N;
                                        }
                                        (*Abuf_row++) += sum;
                                        
                                        WY_row -= k_ub;
                                        A_row -= A_diff_val;
                                    }
                                    Abuf_row += Abuf_dif_value;
                                    WY_dstripe_shifted += N;
                                    A_row -= j_ub;
                                }
                            }
                        }
                    }
                }
            }
        }

        double* A_shifted_toUpdate = A + d1_shift*N + d2_shift;
        const double* Abuf_shifted_toUpdate = Abuf + (A_shifted_toUpdate - A);
        for (int i = d1_shift; i < N; ++i)
        {
            memcpy(A_shifted_toUpdate, Abuf_shifted_toUpdate, copy_length);
            A_shifted_toUpdate += N;
            Abuf_shifted_toUpdate += N;
        }
    }// WHILE

    free(work_memory);

    return A;
}

double* QR_WY_standard(double* A, const int N, const int r)
{
    double* work_memory = reinterpret_cast<double*>(
        calloc(N + r + 2*N*r + 2*N*N, sizeof(double)));

    double* v    = new(work_memory) double[N];
    double* w    = new(v + N)       double[r];
    double* W    = new(w + r)       double[N * r];
    double* Y    = new(W + N * r)   double[N * r];
    double* WY   = new(Y + N * r)   double[N * N];
    double* Abuf = new(WY + N * N)  double[N * N];

    // lambda указывает на начало текущего блока
    for (int lambda = 0; lambda < N; lambda += r)
    {
        // Указывает на начало следующего блока за текущим
        const int t = min(lambda + r, N);

        // Выполнение преобразований над текущим блоком
        for (int it = lambda; it < t; ++it)
        {
            double norm = 0.0;
            double scalar = 1.0;
            double buf;
            // * Вычисление вектора Хаусхолдера
            for (int j = it; j < N; ++j)
            {
                buf = A[j*N + it];
                v[j] = buf;
                norm += buf * buf;
            }
            const double A_diag_el = A[it * N + it];
            const double diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + diag_el_sign * sqrt(norm);
            // используется нормировка, т.ч. v[it] = 1
            for (int j = it + 1; j < N; ++j)
            {
                buf = (v[j] /= beta);
                scalar += buf * buf;
            }
            v[it] = 1.0;

            // * Применение преобразования
            beta = -2.0 / scalar;
            // Вычисление вектора w = beta * (A(:,lambda:t-1)^t) * v .
            // начинаем с it, т.к. it-ая матрица Хаусхолдера действует на матрицу A(it:,it:)
            for (int j = it; j < t; ++j)    
            {
                double sum = 0.0;
                for (int k = it; k < N; ++k)
                {
                    sum += A[k*N + j] * v[k];
                }
                w[j-it] = beta * sum;
            }
            // Преобразование текущего блока матрицы A
            for (int j = it; j < N; ++j)
            {
                double* Aj = A + j*N;
                const double vj = v[j];
                for (int k = it; k < t; ++k)
                {
                    Aj[k] += vj * w[k-it];
                }
            }

            // * Заполнение поддиагнальной части it-го столбца матрицы A
            //   информативной частью вектора Хаусхолдера
            for (int j = it + 1; j < N; ++j)
            {
                A[j*N + it] = v[j];
            }
        }// FOR
        
        // * Вычисление WY-представления произведения матриц Хаусхолдера
        const int d = t - lambda;
        memset(Y + lambda*r, 0, r * min(N - lambda, r) * sizeof(double));
        for (int it = 0; it < d; ++it)
        {
            const int shift = lambda + it;
            double scalar = 1.0, buf;
            memset(v + lambda, 0, it*sizeof(double));
            // * Вычисление нормы и ск. произведения вектора Хаусхолдера
            for (int j = shift + 1; j < N; ++j)
            {
                buf = A[j*N + shift];
                scalar += buf * buf;
                v[j] = buf;
            }
            v[shift] = 1.0;
            double beta = -2.0 / scalar;

            // it - число накопленных столбцов в матрицах W и Y (= номеру дозаписываемого столбца)
            if (it == 0)
            {
                // Начальное заполнение первого столбца
                // замечание: первые lambda строк матриц W и Y - нулевые,
                //            матрицы имеют ступенчатый вид.
                for (int i = lambda; i < N; ++i)
                {
                    Y[i*r] = v[i];
                    W[i*r] = beta * v[i];
                }
            }
            else
            {
                // Вычисление произведения (Y^t) * v
                for (int i = 0; i < it; ++i)
                {
                    double sum = 0.0;
                    for (int j = shift; j < N; ++j)
                    {
                        sum += Y[j*r + i] * v[j];
                    }
                    w[i] = sum;
                }

                // Вычисление произведения W * ((Y^t) * v)
                for (int i = lambda; i < N; ++i)
                {
                    double* Wi = W + i*r;
                    double sum = 0.0;
                    for (int j = 0; j < it; ++j)
                    {
                        sum += Wi[j] * w[j];
                    }

                    buf = v[i];
                    Wi[it] = beta * (buf + sum);
                    Y[Wi + it - W] = buf;
                }
            }
        }

        // * Преобразование остальной части матрицы A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)
        
        // Умножение Y * W^t
        for (int i = lambda; i < N; ++i)
        {
            double* WYi = WY + i*N;
            const double* Yi = Y + i*r;
            for (int j = lambda; j < N; ++j)
            {
                const double* Wj = W + j*r;
                double sum = 0.0;
                for (int k = 0; k < d; ++k)
                {
                    sum += Yi[k] * Wj[k];
                }
                WYi[j] = sum;
            }
        }

        // Умножение A(lamda:M-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1)
        for (int i = lambda; i < N; ++i)
        {
            double* Abuf_i = Abuf + i*N;
            const double* Ai = A + i*N;
            const double* WYi = WY + i*N;
            for (int j = t; j < N; ++j)
            {
                double sum = Ai[j];
                for (int k = lambda; k < N; ++k)
                {
                    sum += WYi[k] * A[k*N+j];
                }
                Abuf_i[j] = sum;
            }
        }

        // копирование преобразованной части матрицы A
        for (int i = lambda; i < N; ++i)
        {
            memcpy(A + i*N + t,
                   Abuf + i*N + t,
                   (N - t) * sizeof(double));
        }
    }// WHILE

    free(work_memory);

    return A;
}

double* QR_WY_block(double* A, const TaskData& parameters)
{
    const int N = parameters.M_ROWS;
    const int b1 = parameters.B_ROWS;
    const int b2 = parameters.B_COLS;
    const int block_count_in_col = static_cast<int>(ceil(1.0 * N / b1));
    const int block_count_in_row = static_cast<int>(ceil(1.0 * N / b2));

    double* work_memory = reinterpret_cast<double*>(
        calloc(N + b2 + 2*N*b2 + 2*N*N, sizeof(double)));

    double* v =       new(work_memory)  double[N];
    double* w =       new(v + N)        double[b2];
    double* W =       new(w + b2) double[N*b2];
    double* Y =       new(W + N*b2)     double[N*b2];
    double* WY =      new(Y + N*b2)     double[N*N];
    double* Abuf =    new(WY + N*N)     double[N*N];

    int curr_i_block_ind = -1;
    int curr_j_block_ind = -1;
    // lambda указывает на начало текущего блока
    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // координаты блока, в котором находится элемент (lambda,lambda).
        // замечание: i-координата может измениться из-за разных размеров b1 и b2 блоков. 
        curr_i_block_ind = lambda / b1;
        ++curr_j_block_ind;

        // Указывает на начало следующего по горизонтали блока
        const int t = min(lambda + b2, N);
        // оставшееся число блоков справа от текущего по строке
        const int row_b_count = block_count_in_row - curr_j_block_ind - 1;
        // оставшееся число блоков включая текущий и ниже по столбцу
        const int col_b_count = block_count_in_col - curr_i_block_ind;
        // сдвиги по вертикали и горизонтали до текущего блока
        const int d1_shift = curr_i_block_ind * b1;
        const int d2_shift = curr_j_block_ind * b2;

        // действительная величина ширины текущего блочного столбца
        const int d2 = t - lambda;
        
        // * Выполнение преобразований над текущим блочным столбцом
        for (int it = 0; it < d2; ++it)
        {
            // сдвиг до 'lambda+it'-го элемнта в столбце
            const int full_shift = lambda + it;
            // фиксация перехода в следующий по вертикали блок
            const int b_begin = (full_shift - d1_shift) / b1;

            memset(w, 0, d2 * sizeof(double));
            memset(v + lambda, 0, it * sizeof(double));
            
            // * Вычисление вектора Хаусхолдера
            double norm(0.0);
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                // сдвиг линии блока от линии матрицы
                const int mb_shift = d1_shift + jb * b1;
                // верхняя и нижняя границы чтения строк в блоке jb (локальная нумерация)
                const int j_lb = max(full_shift, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // высота блока jb
                const int block_h = min(b1, N - mb_shift);
                // указатель на 'it'-ый элемент 'j_lb'-ой строки блока 'jb'
                const double* Ablock_shifted = A + mb_shift*N +
                    block_h*d2_shift + j_lb*d2 + it;
                double* v_shifted = v + mb_shift + j_lb;
                // j - номер строки в блоке jb (локальная нумерация)
                for (int j = 0; j < j_ub - j_lb; ++j)
                {
                    double buf(*Ablock_shifted);
                    v_shifted[j] = buf;
                    norm += buf * buf;
                    Ablock_shifted += d2;
                }
            }

            const int diag_el_block_h = min(b1, N - (d1_shift + b_begin*b1));
            const double A_diag_el = *(A + d1_shift*N + b_begin*b1*N +
                d2_shift*diag_el_block_h + (full_shift % b1)*d2 + it);
            const double diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + diag_el_sign * sqrt(norm);
            // используется нормировка, т.ч. v[it] = 1
            double scalar(1.0);
            for (int j = full_shift + 1; j < N; ++j)
            {
                double buf(v[j] /= beta);
                scalar += buf * buf;
            }
            v[full_shift] = 1.0;
            beta = -2.0 / scalar;
        
            // * Применение преобразования
            // Вычисление вектора w = beta * (A(:,lambda:t-1)^t) * v .
            // Начинаем с it, т.к. it-ая матрица Хаусхолдера действует на матрицу A(it:,it:)
            for (int kb = b_begin; kb < col_b_count; ++kb)
            {
                // сдвиг линии блока от линии матрицы
                const int mb_shift = d1_shift + kb * b1;
                // верхняя и нижняя границы чтения строк в блоке jb (локальная нумерация)
                const int k_lb = max(full_shift, mb_shift) - mb_shift;
                const int k_ub = min(mb_shift + b1, N) - mb_shift;
                const int k_iter_count = k_ub - k_lb;
                // Величины сдвигов указателей на 'A' и на 'v'
                // с возвратом их на исходные позиции
                // для выполнения очереной итерации умножения
                const int A_diff_val = k_iter_count * d2 - 1;
                // высота блока jb
                const int block_h = min(b1, N - mb_shift);
                // указатель на начало блока jb
                const double* Ablock_shifted = A + mb_shift*N +
                    block_h*d2_shift + k_lb*d2 + it;
                const double* v_shifted = v + mb_shift + k_lb;
                // j - номер столбца в блоке kb (локальная нумерация)
                for (int j = it; j < d2; ++j)
                {
                    double sum(0.0);
                    for (int k = 0; k < k_iter_count; ++k)
                    {
                        sum += (*Ablock_shifted) * v_shifted[k];
                        Ablock_shifted += d2;
                    }
                    w[j] += sum;

                    Ablock_shifted -= A_diff_val;
                }
            }
            for (int k = it; k < d2; ++k)
            {
                w[k] *= beta;
            }

            // Преобразование текущего блока матрицы A
            const double* w_shifted_TRANSF = w + it;
            const int w_diff_val_TRANSF = d2 - it;
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                // сдвиг линии блока от линии матрицы
                const int mb_shift = d1_shift + jb*b1;
                // верхняя и нижняя границы чтения строк в блоке jb (локальная нумерация)
                const int j_lb = max(full_shift, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // высота блока jb
                const int j_block_h = min(b1, N - mb_shift);
                // указатель на 'it'-ый элемент 'j_lb'-ой строки блока 'jb'
                double* Arow_shifted = A + mb_shift*N +
                    j_block_h*d2_shift + j_lb*d2 + it;
                const double* v_shifted = v + mb_shift + j_lb;
                // j - номер строки в блоке jb (локальная нумерация)
                for (int j = 0; j < j_ub - j_lb; ++j)
                {
                    const double vj(v_shifted[j]);
                    // k - номер столбца в блоке jb (локальная нумерация)
                    for (int k = 0; k < w_diff_val_TRANSF; ++k)
                    {
                        (*Arow_shifted++) += vj * (*w_shifted_TRANSF++);
                    }
                    Arow_shifted += it;
                    w_shifted_TRANSF -= w_diff_val_TRANSF;
                }
            }
            
            // * Заполнение поддиагнальной части 'it'-го столбца матрицы 'A'
            //   информативной частью вектора Хаусхолдера
            const int b_inc_begin = (full_shift+1 - curr_i_block_ind*b1) / b1;
            for (int jb = b_inc_begin; jb < col_b_count; ++jb)
            {
                // сдвиг линии блока от линии матрицы
                const int mb_shift = d1_shift + jb*b1;
                // верхняя и нижняя границы чтения строк в блоке jb (локальная нумерация)
                const int j_lb = max(full_shift+1, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // высота блока jb
                const int j_block_h = min(b1, N - mb_shift);
                // указатель на начало блока jb
                double* Ablock_shifted = A + mb_shift*N +
                    j_block_h*d2_shift + j_lb*d2 + it;
                const double* v_shifted = v + mb_shift + j_lb;
                // j - номер строки в блоке jb (локальная нумерация)
                for (int j = 0; j < j_ub - j_lb; ++j)
                {
                    *Ablock_shifted = v_shifted[j];
                    Ablock_shifted += d2;
                }
            }
        }// FOR

        memset(W + lambda*b2, 0, (N - lambda)*b2*sizeof(double));
        memset(Y + lambda*b2, 0, (N - lambda)*b2*sizeof(double));

        // * Вычисление WY-представления произведения матриц Хаусхолдера
        for (int it = 0; it < d2; ++it)
        {
            // сдвиг до 'lambda+it'-го элемнта в столбце
            const int full_shift = lambda + it;
            // фиксация перехода в следующий по вертикали блок
            const int b_begin = (full_shift - d1_shift) / b1;

            memset(w, 0, d2*sizeof(double));
            memset(v + lambda, 0, it*sizeof(double));

            double scalar(1.0), buf(0.0);
            // * Вычисление нормы и ск. произведения вектора Хаусхолдера
            for (int jb = (full_shift+1 - d1_shift)/b1; jb < col_b_count; ++jb)
            {
                // сдвиг линии блока от линии матрицы
                const int mb_shift = d1_shift + jb * b1;
                // верхняя и нижняя границы чтения строк в блоке jb (локальная нумерация)
                const int j_lb = max(full_shift+1, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // высота блока jb
                const int block_h = min(b1, N - mb_shift);
                // указатель на начало блока jb
                const double* Arow_shifted = A + mb_shift*N +
                    block_h*d2_shift + j_lb*d2 + it;
                double* v_shifted = v + mb_shift + j_lb;
                // j - номер строки в блоке jb (локальная нумерация)
                for (int j = 0; j < j_ub - j_lb; ++j)
                {
                    buf = *Arow_shifted;
                    scalar += buf * buf;
                    v_shifted[j] = buf;
                    Arow_shifted += d2;
                }
            }
            v[full_shift] = 1.0;
            double beta = -2.0 / scalar;
            
            // 'it' - число накопленных столбцов в матрицах 'W' и 'Y'
            // (= номеру дописываемого столбца)
            if (it == 0)
            {
                // Начальное заполнение первого столбца
                // замечание: первые 'lambda' строк матриц 'W' и 'Y' - нулевые,
                //            матрица 'Y' имеет ступенчатый вид.
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    // сдвиг линии блока от линии матрицы
                    const int mb_shift = d1_shift + ib * b1;
                    // верхняя и нижняя границы чтения строк в блоке ib (локальная нумерация)
                    const int i_lb = max(full_shift, mb_shift) - mb_shift;
                    const int i_ub = min(mb_shift + b1, N) - mb_shift;
                    // указатели на начала блоков ib матриц Y и W
                    double* Yrow = Y + (mb_shift + i_lb)*b2;
                    double* Wrow = W + (Yrow - Y);
                    const double* v_shifted = v + mb_shift + i_lb;
                    for (int i = 0; i < i_ub - i_lb; ++i)
                    {
                        double v_buf(v_shifted[i]);
                        *Yrow = v_buf;
                        *Wrow = beta * v_buf;
                        Yrow += b2;
                        Wrow += b2;
                    }
                }
            }
            else
            {
                // Вычисление произведения (Y^t) * v
                double* wp_Yv = w;
                const int W_Y_diff_val = b2 - it;
                for (int jb = b_begin; jb < col_b_count; ++jb)
                {
                    // сдвиг линии блока от линии матрицы
                    const int mb_shift = d1_shift + jb * b1;
                    // верхняя и нижняя границы чтения строк в блоке jb (глобальная нумерация)
                    const int j_lb = max(full_shift, mb_shift);
                    const int j_ub = min(mb_shift + b1, N);
                    // указатель на начало блока jb
                    const double* Yrow = Y + j_lb*b2;
                    double* v_shifted = v + j_lb;
                    // j - номер строки в блоке jb (глобальная нумерация)
                    for (int j = j_lb; j < j_ub; ++j)
                    {
                        const double vj = *v_shifted++;
                        for (int i = 0; i < it; ++i)
                        {
                            (*wp_Yv++) += (*Yrow++) * vj;
                        }
                        Yrow += W_Y_diff_val;
                        wp_Yv -= it;
                    }
                }
                
                // Вычисление произведения W * ((Y^t) * v) <==> W * w
                const double* wp_Ww = w;
                for (int ib = 0; ib < col_b_count; ++ib)
                {
                    // сдвиг линии блока от линии матрицы
                    const int mb_shift = d1_shift + ib * b1;
                    // верхняя и нижняя границы чтения строк в блоке ib (локальная нумерация)
                    const int i_lb = max(lambda, mb_shift);
                    const int i_ub = min(mb_shift + b1, N);
                    // указатели на начала блоков ib матриц Y и W
                    double* Yrow_shifted = Y + i_lb*b2 + it;
                    double* Wrow = W + i_lb*b2;
                    const double* v_shifted = v + i_lb;
                    for (int i = 0; i < i_ub - i_lb; ++i)
                    {
                        double sum(0.0);
                        for (int j = 0; j < it; ++j)
                        {
                            sum += Wrow[j] * wp_Ww[j];
                        }
                        buf = v_shifted[i];
                        Wrow[it] = beta * (buf + sum);
                        *Yrow_shifted = buf;

                        Wrow += b2;
                        Yrow_shifted += b2;
                    }
                }
            }
        }

        // * Преобразование остальной части матрицы A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        // Умножение Y * W^t
        // Матрица 'WY' будет иметь 'b1' x 'b1' размещение
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib * b1;
            const int i_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int i_ub = min(b1, N - i_mb_shift);
            const int i_iter_count = i_ub - i_lb;
            const int Y_diff_val = (i_ub - i_lb) * b2;
            const int ib_stripe_h = min(b1, N - i_mb_shift);
            const double* Yrow = Y + (i_mb_shift + i_lb) * b2;
            double* WYstripe = WY + i_mb_shift*N + ib_stripe_h*d1_shift;

            for (int jb = 0; jb < col_b_count; ++jb)
            {
                const int j_mb_shift = d1_shift + jb * b1;
                const int j_lb = max(lambda, j_mb_shift) - j_mb_shift;
                const int j_ub = min(b1, N - j_mb_shift);
                const int j_iter_count = j_ub - j_lb;
                const int jb_block_w = min(b1, N - j_mb_shift);
                const int WY_diff_val = jb_block_w - (j_ub - j_lb);
                const int W_diff_val = (j_ub - j_lb) * b2;
                const double* Wrow = W + (j_mb_shift + j_lb) * b2;
                double* WYrow_shifted = WYstripe + jb * ib_stripe_h * b1 +
                    i_lb * jb_block_w + j_lb;

                for (int i = 0; i < i_iter_count; ++i)
                {
                    // Ограничение, связанное со ступенчатым видом матрицы Y
                    const int kbound = min(t, i_mb_shift+i_lb + i+1) - lambda;
                    const int Wrow_diff_val = b2 - kbound;
                                        
                    for (int j = 0; j < j_iter_count; ++j)
                    {
                        double sum(0.0);
                        for (int k = 0; k < kbound; ++k)
                        {
                            sum += (*Yrow++) * (*Wrow++);
                        }
                        (*WYrow_shifted++) = sum;

                        Yrow -= kbound;
                        Wrow += Wrow_diff_val;
                    }
                    WYrow_shifted += WY_diff_val;
                    Yrow += b2;
                    Wrow -= W_diff_val;
                }
                Yrow -= Y_diff_val;
            }
        }

        copy_minor_block(Abuf, A, parameters, d1_shift, d2_shift);

        // Умножение A(lamda:M-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1)
        // Размещение: (b1 x b2) = (b1 x b1) * (b1 x b2)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib * b1;
            const int i_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int ib_stripe_h = min(b1, N - i_mb_shift);
            const int i_iter_count = ib_stripe_h - i_lb;
            const double* A_stripe_shifted = A + i_mb_shift*N + t*ib_stripe_h;
            double* Abuf_stripe_shifted = Abuf + (A_stripe_shifted - A);
            const double* WY_stripe_shifted = WY + i_mb_shift*N + d1_shift*ib_stripe_h;

            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_mb_shift = t + jb * b2;
                const int jb_block_w = min(b2, N - j_mb_shift);
                const int Abuf_diff_val = (ib_stripe_h - i_lb) * jb_block_w;
                double* Abuf_row = Abuf_stripe_shifted +
                    jb * ib_stripe_h * b2 + i_lb * jb_block_w;

                for (int kb = 0; kb < col_b_count; ++kb)
                {
                    const int k_mb_shift = d1_shift + kb * b1;
                    const int k_lb = max(lambda, k_mb_shift) - k_mb_shift;
                    const int kb_block_size = min(b1, N - k_mb_shift);
                    const int k_iter_count = kb_block_size - k_lb;
                    const int Arow_diff_val = k_iter_count * jb_block_w - 1;
                    const double* WYrow_shifted = WY_stripe_shifted +
                        kb * ib_stripe_h * b1 + i_lb * kb_block_size + k_lb;
                    const double* Arow_shifted = A + k_mb_shift * N +
                        j_mb_shift * kb_block_size + k_lb * jb_block_w;

                    for (int i = 0; i < i_iter_count; ++i)
                    {
                        for (int j = 0; j < jb_block_w; ++j)
                        {
                            double sum(0.0);
                            for (int k = 0; k < k_iter_count; ++k)
                            {
                                sum += WYrow_shifted[k] * (*Arow_shifted);
                                Arow_shifted += jb_block_w;
                            }
                            (*Abuf_row++) += sum;

                            Arow_shifted -= Arow_diff_val;
                        }
                        WYrow_shifted += kb_block_size;
                        Arow_shifted -= jb_block_w;
                    }
                    Abuf_row -= Abuf_diff_val;
                }
            }
        }

        // Копирование преобразованной части матрицы A
        copy_minor_block(A, Abuf, parameters, d1_shift, d2_shift);

    }// WHILE

    free(work_memory);

    return A;
}

double* QR_WY_double_block(double* A, const TaskData& parameters)
{
    // Main parameters of data layout
    const int N = parameters.M_ROWS;
    const int b1 = parameters.B_ROWS;
    const int b2 = parameters.B_COLS;
    const int db1 = parameters.D_ROWS;
    const int db2 = parameters.D_COLS;

    const int bsizes_ratio = static_cast<int>(ceil(1.0 * b2 / b1));
    // Big blocks counts in matrix
    const int block_count_in_col = static_cast<int>(ceil(1.0 * N / b1));
    const int block_count_in_row = static_cast<int>(ceil(1.0 * N / b2));
    // Small blocks counts in big blocks
    const int dblock_count_in_bcol = static_cast<int>(ceil(1.0 * b1 / db1));
    const int dblock_count_in_brow = static_cast<int>(ceil(1.0 * b2 / db2));
    const int dblock_count_in_diff_bcol = (N % b1 == 0) ?
        dblock_count_in_bcol :
        static_cast<int>(ceil(1.0 * (N % b1) / db1));
    const int dblock_count_in_diff_brow = (N % b2 == 0) ?
        dblock_count_in_brow :
        static_cast<int>(ceil(1.0 * (N % b2) / db2));

    double* work_memory = reinterpret_cast<double*>(
            calloc(N + b2 + 2*N*b2 + 2*N*N + db1, sizeof(double)));

    double* v =       new(work_memory)   double[N];
    double* w =       new(v + N)         double[b2];
    double* sum_vec = new(w + b2)        double[db1];
    double* W =       new(sum_vec + db1) double[N*b2];
    double* Y =       new(W + N*b2)      double[N*b2];
    double* WY =      new(Y + N*b2)      double[N*N];
    double* Abuf =    new(WY + N*N)      double[N*N];

    int curr_i_block_ind = -1;
    int curr_j_block_ind = -1;

    // Layout data for WY-matrix
    TaskClass WY_task_class;

    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // Indices of big block, that contains element (lambda,lambda).
        // Remark: row-index can be not constant when b1 != b2.
        curr_i_block_ind = lambda / b1;
        ++curr_j_block_ind;

        // Index of first column of next block-column
        const int t = min(lambda + b2, N);
        // Count of big blocks from the right of current big block in its block-stripe
        const int row_b_count = block_count_in_row - curr_j_block_ind - 1;
        // Count of big blocks below and including current big block in its block-column
        const int col_b_count = block_count_in_col - curr_i_block_ind;
        // Horizontal and vertical shifts to diagonal block
        const int d1_shift = curr_i_block_ind * b1;
        const int d2_shift = curr_j_block_ind * b2;
        // Real value of width of current block-column
        const int d2 = t - lambda;
        // Count of small blocks in current big block in row-direction
        const int dblock_count_in_curr_block_row = 
            (curr_j_block_ind == block_count_in_row - 1) ?
            dblock_count_in_diff_brow : dblock_count_in_brow;
        
        // * Execution of current block-column transformations...
        for (int it = 0; it < d2; ++it)
        {
            memset(w, 0, d2 * sizeof(double));
            memset(v + lambda, 0, it * sizeof(double));

            // Shift for 'lambda+it'-th element of matrix column
            const int full_shift = lambda + it;
            // Индекс большого блока, содержащего 'lambda+it'-ый диагональный элемент,
            // считая от блока, содержащего 'lambda'-ый элемент
            const int b_begin = (full_shift - d1_shift) / b1;
            // Index of small-block-column,
            // which contains 'lambda+it'-th column
            const int dblock_col_index = it / db2;
            // Index of small-block-row of diagonal big block,
            // which contains 'full_shift'-th row of matrix
            const int dblock_row_index = (full_shift - (d1_shift + b_begin*b1)) / db1;
            // Width of this small block
            const int dblock_width = min(db2, d2 - dblock_col_index*db2);
            const int loc_it = it - dblock_col_index*db2;

            // * Householder vector computation...
            // Using of 'it'-th column of current block-column,
            // starting with (lambda+it)-th position
            double norm(0.0);
            for (int ib = b_begin; ib < col_b_count; ++ib)
            {
                // Shift of current big block first line
                // from matrix first line
                const int mb_shift = d1_shift + ib * b1;
                // Upper and lower bounds for 'ib'-th block lines reading.
                // 'i_lb' is required, because we need to use block rows
                // starting with 'full_shift'-th matrix row,
                // that can be mapped not only to first block row.
                const int i_lb = max(full_shift, mb_shift) - mb_shift;  // local indexation
                // 'i_ub' is required, because last matrix row,
                // that can be mapped not only to last block row, but smaller.
                const int i_ub = min(mb_shift + b1, N);  // global indexation
                // 'ib'-th block height (number of rows)
                const int ib_block_h = min(b1, N - mb_shift);
                // 'ib'-th block beginning pointer
                const double* Ablock = A + mb_shift * N + ib_block_h * d2_shift;

                // Upper and lower bounds for small blocks stripes reading.
                // Reasons of using these bounds are the same as 'i_lb' and 'i_ub'.
                const int id_lb = (i_lb == 0) ? 0 : (i_lb / db1);
                const int id_ub = (ib == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int loc_shift = id * db1;
                    // Shift of current small-block-row from matrix first row
                    const int bd_shift = mb_shift + loc_shift;
                    // Bounds for current small block reading.
                    // Reasons of using are the same as 'i_lb' and 'i_ub'
                    const int lb = max(full_shift, bd_shift) - bd_shift;
                    const int ub = min(bd_shift + db1, i_ub) - bd_shift;
                    // 'id'-th small block height (number of rows)
                    const int id_block_h = min(db1, ib_block_h - loc_shift);
                    // Pointer to beginning of required small block in 'id'-th stripe
                    const double* Adblock_shifted = Ablock + loc_shift*d2 +
                        dblock_col_index*id_block_h*db2 + lb*dblock_width + loc_it;
                    double* v_shifted = v + bd_shift + lb;
                    
                    for (int i = 0; i < ub - lb; ++i)
                    {
                        double buf(*Adblock_shifted);
                        v_shifted[i] = buf;
                        norm += buf * buf;

                        Adblock_shifted += dblock_width;
                    }
                }
            }
                        
            // Diagonal big block height
            const int diag_el_block_h = min(b1, N - (d1_shift + b_begin*b1));
            // Diagonal small block height
            const int diag_el_dblock_h = min(db1, diag_el_block_h - dblock_row_index*db1);
            const double A_diag_el = *(A + d1_shift*N + b_begin*b1*N +  // shift for big block stripe
                d2_shift * diag_el_block_h +  // shift for big block beginning in its stripe
                ((full_shift % b1) / db1)*db1*d2 +  // shift for small block stripe
                (dblock_col_index*diag_el_dblock_h*db2) +   // shift for small block beginning in its stripe
                ((full_shift % b1) % db1)*dblock_width + loc_it);  // shift for diagonal element in small block

            const double diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + diag_el_sign * sqrt(norm);
            // Use normalization such that v[it] = 1
            double scalar(1.0);
            double* v_shifted_NORM = v + full_shift + 1;
            for (int j = 0; j < N - full_shift - 1; ++j)
            {
                double buf(v_shifted_NORM[j] /= beta);
                scalar += buf * buf;
            }
            v[full_shift] = 1.0;
            beta = -2.0 / scalar;

            // * Transformation applying...

            // Computation of vector w = beta * (A(:,lambda:t-1)^t) * v .
            // Starting with 'full_shift'-th row and column of matrix,
            // because 'full_shift'-th Householder matrix really operates
            // with A(full_shift:*, full_shift:t-1).
            // Pattern of data access is as follows:
            // 1. Fix big block (kb)
            // 2. Fix small block stripe in fixed big block (kd)
            // 3. Fix small block in fixed stripe (jd)
            // 4. Fix column in fixed small block (j)
            // 5. Produce passing of fixed column (k)
            for (int kb = b_begin; kb < col_b_count; ++kb)
            {
                const int mb_shift = kb*b1 + d1_shift;
                const int k_lb = max(full_shift, mb_shift) - mb_shift;
                const int k_ub = min(mb_shift + b1, N);
                const int block_h = min(b1, N - mb_shift);
                const double* Ablock = A + mb_shift*N + block_h*d2_shift;
                const double* v_kb_shift = v + mb_shift;
                
                const int kd_lb = (k_lb == 0) ? 0 : k_lb / db1;
                const int kd_ub = (kb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int kd = kd_lb; kd < kd_ub; ++kd)
                {
                    const int loc_shift = kd * db1;
                    const int bd_shift = loc_shift + mb_shift;
                    const int lb = max(full_shift, bd_shift) - bd_shift;
                    const int ub = min(db1, k_ub - bd_shift);
                    const int k_iter_count = ub - lb;
                    const int dblock_h = min(db1, block_h - loc_shift);
                    // Pointer to beginning of 'kb'-th small block stripe
                    const double* Adstripe = Ablock + loc_shift*d2;
                    const double* v_kd_shifted = v_kb_shift + loc_shift + lb;
                    
                    // const int& jd_lb = dblock_col_index;
                    // const int& jd_ub = dblock_count_in_curr_block_row;
                    for (int jd = dblock_col_index; jd < dblock_count_in_curr_block_row; ++jd)
                    {
                        const int loc_drow_shift = jd * db2;
                        // Width of this small block
                        const int dblock_w = min(db2, d2 - loc_drow_shift);
                        const int j_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int j_ub = min(d2 - loc_drow_shift, db2);
                        const int j_iter_count = j_ub - j_lb;
                        const int A_diff_val = k_iter_count * dblock_w - 1;
                        // Pointer to beginning of 'jb'-th small block
                        const double* Arow_shifted = Adstripe +
                            dblock_h * loc_drow_shift + lb * dblock_w + j_lb;
                        double* w_shifted = w + loc_drow_shift + j_lb;

                        for (int j = 0; j < j_iter_count; ++j)
                        {
                            double sum(0.0);
                            for (int k = 0; k < k_iter_count; ++k)
                            {
                                sum += (*Arow_shifted) * v_kd_shifted[k];
                                Arow_shifted += dblock_w;
                            }
                            Arow_shifted -= A_diff_val;
                            w_shifted[j] += sum;
                        }
                    }
                }
            }
            for (int k = it; k < d2; ++k)
            {
                w[k] *= beta;
            }

            // Transformation of current block column of matrix A:
            // A(fullshift:N, it:t-1) += v[fullshift:N] * (w[it:t-1]^t)
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                const int mb_shift = d1_shift + jb * b1;
                const int jb_lb = max(full_shift, mb_shift) - mb_shift;
                const int jb_ub = min(mb_shift + b1, N);
                const int jb_block_h = min(b1, N - mb_shift);
                double* Ablock = A + mb_shift*N + jb_block_h*d2_shift;

                const int jd_lb = (jb_lb == 0) ? 0 : jb_lb / db1;
                const int jd_ub = (jb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int jd = jd_lb; jd < jd_ub; ++jd)
                {
                    const int loc_shift = jd * db1;
                    const int bd_shift = loc_shift + mb_shift;
                    const int j_lb = max(full_shift, bd_shift) - bd_shift;
                    const int j_ub = min(db1, jb_ub - bd_shift);
                    const int j_it_count = j_ub - j_lb;
                    const int jd_block_h = min(db1, jb_block_h - loc_shift);
                    double* Adstripe = Ablock + loc_shift*d2;
                    const double* v_shifted = v + bd_shift + j_lb;

                    for (int kd = dblock_col_index; kd < dblock_count_in_curr_block_row; ++kd)
                    {
                        const int loc_drow_shift = kd * db2;
                        const int k_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int dblock_w = min(db2, d2 - loc_drow_shift);
                        const int k_iter_count = dblock_w - k_lb;
                        double* Arow_shifted = Adstripe + jd_block_h * loc_drow_shift +
                            j_lb * dblock_w + k_lb;
                        const double* w_shifted = w + loc_drow_shift + k_lb;

                        for (int j = 0; j < j_it_count; ++j)
                        {
                            const double vj = v_shifted[j];
                            for (int k = 0; k < k_iter_count; ++k)
                            {
                                (*Arow_shifted++) += vj * (*w_shifted++);
                            }

                            w_shifted -= k_iter_count;
                            Arow_shifted += k_lb;
                        }
                    }
                }
            }

            // Filling up of underdiagonal part of matrix 'full_shift'-th column
            // with informative part of Householder vector 'v'
            const int b_inc_begin = (full_shift + 1 - d1_shift) / b1;
            for (int ib = b_inc_begin; ib < col_b_count; ++ib)
            {
                const int mb_shift = d1_shift + ib*b1;
                const int ib_lb = max(full_shift + 1, mb_shift) - mb_shift;
                const int ib_ub = min(mb_shift + b1, N);
                const int ib_block_h = min(b1, N - mb_shift);
                double* Ablock = A + mb_shift*N + ib_block_h*d2_shift;

                const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
                const int id_ub = (ib == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int loc_shift = id * db1;
                    const int bd_shift = loc_shift + mb_shift;
                    const int i_lb = max(full_shift + 1, bd_shift) - bd_shift;
                    const int i_ub = min(db1, ib_ub - bd_shift);
                    const int id_block_h = min(db1, ib_block_h - loc_shift);
                    double* Arow_shifted = Ablock + loc_shift*d2 +
                        dblock_col_index*id_block_h*db2 +
                        i_lb*dblock_width + loc_it;
                    const double* v_shifted = v + bd_shift + i_lb;
                    
                    for (int i = 0; i < i_ub - i_lb; ++i)
                    {
                        *Arow_shifted = v_shifted[i];
                        Arow_shifted += dblock_width;
                    }
                }
            }
        }// FOR

        // * Computation of WY-representation of Householder matrices production...
        memset(W + d1_shift*b2, 0, min(bsizes_ratio*b1, N - d1_shift)*b2*sizeof(double));
        memset(Y + d1_shift*b2, 0, min(bsizes_ratio*b1, N - d1_shift)*b2*sizeof(double));
        for (int it = 0; it < d2; ++it)
        {
            memset(w, 0, d2 * sizeof(double));
            memset(v + lambda, 0, it * sizeof(double));

            const int full_shift = lambda + it;
            const int b_begin = (full_shift - d1_shift) / b1;
            const int dblock_col_index = it / db2;
            const int dblock_width = min(db2, d2 - dblock_col_index * db2);
            // Separated value is required because W and Y matrices are
            // always storing as N x b2 with d1 x d2 small blocks. 
            // So, if actual width of column is less then 'b2',
            // then form of small blocks won't change.
            const int dblock_width_WY = min(db2, b2 - dblock_col_index*db2);
            const int loc_it = it - dblock_col_index*db2;

            double scalar = 1.0;
            const int b_inc_begin = (full_shift - d1_shift) / b1;
            for (int ib = b_inc_begin; ib < col_b_count; ++ib)
            {
                const int mb_shift = d1_shift + ib * b1;
                const int i_lb = max(full_shift + 1, mb_shift) - mb_shift;
                const int i_ub = min(mb_shift + b1, N);
                const int ib_block_h = min(b1, N - mb_shift);
                const double* Ablock = A + mb_shift * N +
                    ib_block_h * d2_shift;

                const int id_lb = (i_lb == 0) ? 0 : (i_lb / db1);
                const int id_ub = (ib == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int loc_shift = id * db1;
                    const int bd_shift = mb_shift + loc_shift;
                    const int lb = max(full_shift + 1, bd_shift) - bd_shift;
                    const int ub = min(db1, i_ub - bd_shift);
                    const int id_block_h = min(db1, ib_block_h - loc_shift);
                    const double* Adblock_shifted = Ablock + loc_shift*d2 +
                        dblock_col_index*id_block_h*db2 +
                        loc_it + lb*dblock_width;
                    double* v_shifted = v + bd_shift + lb;

                    for (int i = 0; i < ub - lb; ++i)
                    {
                        double buf(*Adblock_shifted);
                        Adblock_shifted += dblock_width;
                        scalar += buf * buf;
                        v_shifted[i] = buf;
                    }
                }
            }
            v[full_shift] = 1.0;
            double beta = -2.0 / scalar;
                        
            // 'it' is count of accumulated columns in W and Y matrices
            // (equals index of columns, which is added at current iteration)
            if (it == 0)
            {
                // Initial first column filling
                // Remark 1: first 'lambda' of rows in matrices W and Y are filled with zero,
                //           next 'd2' of rows consist lower triangular matrix.
                // Remark 2: W and Y will be allocated with (B1,B2,D1,D2)-layout.
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    const int mb_shift = ib*b1 + d1_shift;
                    const int i_lb = max(lambda, mb_shift) - mb_shift;
                    const int i_ub = min(mb_shift + b1, N);

                    const int id_lb = (i_lb == 0) ? 0 : (i_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;
                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int loc_shift = id * db1;
                        const int bd_shift = loc_shift + mb_shift;
                        const int lb = max(lambda, bd_shift) - bd_shift;
                        const int ub = min(db1, i_ub - bd_shift);
                        // Width of block is 'b2',
                        // because memory was allocated for complete matrices
                        double* Ydblock = Y + bd_shift*b2 + lb*dblock_width_WY;
                        double* Wdblock = W + (Ydblock - Y);
                        const double* v_shifted = v + bd_shift + lb;

                        for (int i = 0; i < ub - lb; ++i)
                        {
                            double v_buf(v_shifted[i]);
                            *Ydblock = v_buf;
                            *Wdblock = beta * v_buf;
                            Ydblock += dblock_width_WY;
                            Wdblock += dblock_width_WY;
                        }
                    }
                }
            }
            else
            {
                // Computation of production w = (Y^t) * v .
                // Pattern of data access is as follows:
                // 1. Fix big block (ib)
                // 2. Fix small block stripe in fixed big block (id)
                // 3. Fix small block in fixed stripe (jd)
                // 4. Fix row in fixed small block (i) and element of vector 'v'
                // 5. Produce passing of fixed row of Y (j) and vector 'w',
                //    using fixed element of vector 'v'
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    const int mb_shift = ib*b1 + d1_shift;
                    const int ib_lb = max(full_shift, mb_shift) - mb_shift;
                    const int ib_ub = min(mb_shift + b1, N);
                    const int ib_block_h = min(b1, N - mb_shift);
                    const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;

                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int loc_shift = id * db1;
                        const int bd_shift = loc_shift + mb_shift;
                        const int i_lb = max(full_shift, bd_shift) - bd_shift;
                        const int i_ub = min(db1, ib_ub - bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int id_block_h = min(db1, ib_block_h - loc_shift);
                        const double* Ydstripe = Y + bd_shift*b2;
                        const double* v_shifted = v + bd_shift + i_lb;

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it - loc_drow_shift, db2);
                            const int dblock_w = min(db2, b2 - loc_drow_shift);
                            const double* Yrow = Ydstripe +
                                id_block_h * loc_drow_shift + i_lb * dblock_w;
                            double* w_shifted = w + loc_drow_shift;

                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                const double vi = v_shifted[i];
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    w_shifted[j] += Yrow[j] * vi;
                                }
                                Yrow += dblock_w;
                            }
                        }
                    }
                }

                // Computation of production W * ((Y^t) * v) <=> W * w .
                // Pattern of data access is as follows:
                // 1. Fix big block (ib)
                // 2. Fix small block stripe in fixed big block (id)
                // 3. Fix small block in fixed stripe (jd)
                // 4. Fix row in fixed small block (i)
                // 5. Produce passing of fixed row of W (j) and vector 'w'
                const int set_size_SUM_VEC = db1 * sizeof(double);
                for (int ib = 0; ib < col_b_count; ++ib)
                {
                    const int mb_shift = ib*b1 + d1_shift;
                    const int ib_lb = max(lambda, mb_shift) - mb_shift;
                    const int ib_ub = min(mb_shift + b1, N);
                    const int ib_block_h = min(b1, N - mb_shift);

                    const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;
                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int loc_shift = id * db1;
                        const int bd_shift = loc_shift + mb_shift;
                        const int i_lb = max(lambda, bd_shift) - bd_shift;
                        const int i_ub = min(db1, ib_ub - bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int id_block_h = min(db1, ib_block_h - loc_shift);
                        double* Wdstripe = W + bd_shift * b2;
                        double* Ydstripe = Y + (Wdstripe - W);
                        // Vector to contain results of production of W's row and vector w
                        memset(sum_vec, 0, set_size_SUM_VEC);

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it, loc_drow_shift + db2) - loc_drow_shift;
                            const int dblock_w = min(db2, b2 - loc_drow_shift);
                            double* Wrow = Wdstripe +
                                id_block_h * loc_drow_shift + i_lb * dblock_w;
                            const double* w_shifted = w + loc_drow_shift;
                            double* sum_vec_shifted = sum_vec + i_lb;

                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                double sum_loc(0.0);
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    sum_loc += Wrow[j] * w_shifted[j];
                                }
                                // Scalar production of W matrix line and vector w
                                sum_vec_shifted[i] += sum_loc;

                                Wrow += dblock_w;
                            }
                        }
                        // Insertion of computed values into 'it'-th column of Y and W
                        const double* v_shifted = v + bd_shift + i_lb;
                        // Pointers to small blocks containing 'it'-th column of Y and W
                        double* Wdblock_add = Wdstripe + dblock_col_index*id_block_h*db2 +
                            i_lb*dblock_width_WY + loc_it;
                        double* Ydblock_add = Ydstripe + (Wdblock_add - Wdstripe);
                        const double* sum_vec_shifted = sum_vec + i_lb;

                        for (int i = 0; i < i_iter_count; ++i)
                        {
                            double buf(v_shifted[i]);
                            *Wdblock_add = beta * (buf + sum_vec_shifted[i]);
                            *Ydblock_add = buf;
                            Wdblock_add += dblock_width_WY;
                            Ydblock_add += dblock_width_WY;
                        }
                    }
                }
            }
        }

        // * Transformation of remaining block-columns of A...
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        // Constructing layout data for WY-matrix
        WY_task_class.makeData(N, N, b1, b1, db1, db1);
        set_minor_double_block(WY, 0, WY_task_class.getDataRef(), d1_shift, d1_shift);

        // Multiplication Y * W^t .
        // Remark: WY matrix get (B1,B1,D1,D1)-layout.
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib * b1;
            const int ib_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int ib_ub = min(i_mb_shift + b1, N);
            const int ib_stripe_h = min(b1, N - i_mb_shift);
            double* WYstripe = WY + i_mb_shift*N;

            const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
            const int id_ub = (ib == col_b_count - 1) ?
                dblock_count_in_diff_bcol :
                dblock_count_in_bcol;
            for (int jb = 0; jb < col_b_count; ++jb)
            {
                const int j_mb_shift = d1_shift + jb *b1;
                const int jb_lb = max(lambda, j_mb_shift) - j_mb_shift;
                const int jb_ub = min(j_mb_shift + b1, N);
                const int jb_stripe_h = min(b1, N - j_mb_shift);
                double* WYblock = WYstripe + ib_stripe_h*j_mb_shift;

                const int jd_lb = (jb_lb == 0) ? 0 : (jb_lb / db1);
                const int jd_ub = (jb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int i_loc_shift = id * db1;
                    const int i_bd_shift = i_mb_shift + i_loc_shift;
                    const int i_lb = max(lambda, i_bd_shift) - i_bd_shift;
                    const int i_ub = min(db1, ib_ub - i_bd_shift);
                    const int i_iter_count = i_ub - i_lb;
                    const int id_block_h = min(db1, ib_stripe_h - i_loc_shift);
                    const double* Ydstripe = Y + i_bd_shift*b2;
                    double* WYdstripe = WYblock + i_loc_shift*jb_stripe_h;
                    
                    // Upper bound for 'kd'-level small-block-cycle is defined by
                    // lower triangular form of matrix Y.
                    // If 'id'-th small-block-stripe is crossing matrix diagonal,
                    // then upper bound for small-block-stripe passing equals
                    // index of small block, that is containing diagonal element
                    // with maximum index within 'id'-th small-block-stripe height limits.
                    const int kd_ud = (i_bd_shift < t) ?
                        ((min(i_bd_shift + db1, ib_ub) - lambda) / db2 + 1):
                        dblock_count_in_curr_block_row;

                    for (int jd = jd_lb; jd < jd_ub; ++jd)
                    {
                        const int j_loc_shift = jd * db1;
                        const int j_bd_shift = j_mb_shift + j_loc_shift;
                        const int j_lb = max(lambda, j_bd_shift) - j_bd_shift;
                        const int j_ub = min(db1, jb_ub - j_bd_shift);
                        const int j_iter_count = j_ub - j_lb;
                        const int jd_block_h = min(db1, jb_stripe_h - j_loc_shift);
                        const int WYrow_diff_val = jd_block_h - j_iter_count;
                        const int WY_diff_val = i_iter_count * jd_block_h;
                        const double* Wdstripe = W + j_bd_shift*b2;
                        double* WYrow_shifted = WYdstripe + id_block_h*j_loc_shift +
                            i_lb * jd_block_h + j_lb;

                        for (int kd = 0; kd < kd_ud; ++kd)
                        {
                            const int k_loc_shift = kd * db2;
                            const int k_ub = min(d2 - k_loc_shift, db2);
                            const int kd_dblock_w = min(db2, b2 - k_loc_shift);
                            const int Wrow_diff_val = kd_dblock_w - k_ub;
                            const int W_diff_val = j_iter_count * kd_dblock_w;
                            const double* Yrow = Ydstripe +
                                id_block_h * k_loc_shift + i_lb * kd_dblock_w;
                            const double* Wrow = Wdstripe +
                                jd_block_h * k_loc_shift + j_lb*kd_dblock_w;

                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                for (int j = 0; j < j_iter_count; ++j)
                                {
                                    double sum(0.0);
                                    for (int k = 0; k < k_ub; ++k)
                                    {
                                        sum += (*Yrow++) * (*Wrow++);
                                    }
                                    (*WYrow_shifted++) += sum;

                                    Yrow -= k_ub;
                                    Wrow += Wrow_diff_val;
                                }

                                WYrow_shifted += WYrow_diff_val;
                                Yrow += kd_dblock_w;
                                Wrow -= W_diff_val;
                            }

                            WYrow_shifted -= WY_diff_val;
                        }
                    }
                }
            }
        }

        copy_minor_double_block(Abuf, A, parameters, d1_shift, d2_shift);

        // Multiplication A(lamda:N-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1) <=>
        //                A(lamda:N-1, t:N-1) += WY * A(lamda:N-1, t:N-1) .
        // Layout pattern: (b1 x b2, d1 x d2) = (b1 x b1, d1 x d1) * (b1 x b2, d1 x d2).
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int ib_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int ib_ub = min(i_mb_shift + b1, N);
            const int ib_stripe_h = min(b1, N - i_mb_shift);
            const double* WY_stripe = WY + i_mb_shift * N;
            double* Abuf_stripe = Abuf + (WY_stripe - WY);

            const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
            const int id_ub = (ib == col_b_count - 1) ?
                dblock_count_in_diff_bcol :
                dblock_count_in_bcol;
            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_mb_shift = t + jb * b2;
                const int jb_ub = min(j_mb_shift + b2, N);
                const int jb_column_w = min(b2, N - j_mb_shift);
                double* Abuf_block = Abuf_stripe + j_mb_shift * ib_stripe_h;

                const int jd_ub = (jb == row_b_count - 1) ?
                    dblock_count_in_diff_brow :
                    dblock_count_in_brow;
                for (int kb = 0; kb < col_b_count; ++kb)
                {
                    const int k_mb_shift = d1_shift + kb * b1;
                    const int kb_lb = max(lambda, k_mb_shift) - k_mb_shift;
                    const int kb_ub = min(k_mb_shift + b1, N);
                    const int kb_stripe_h = min(b1, N - k_mb_shift);
                    const double* WY_block = WY_stripe + ib_stripe_h * k_mb_shift;
                    const double* A_block = A + k_mb_shift * N +
                        j_mb_shift * kb_stripe_h;

                    const int kd_lb = (kb_lb == 0) ? 0 : (kb_lb / db1);
                    const int kd_ub = (kb == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;
                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int i_loc_shift = id * db1;
                        const int i_bd_shift = i_mb_shift + i_loc_shift;
                        const int i_lb = max(lambda, i_bd_shift) - i_bd_shift;
                        const int i_ub = min(db1, ib_ub - i_bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int id_block_h = min(db1, ib_stripe_h - i_loc_shift);
                        const double* WY_dstripe = WY_block + i_loc_shift*kb_stripe_h;
                        double* Abuf_dstripe = Abuf_block + i_loc_shift*jb_column_w;

                        for (int jd = 0; jd < jd_ub; ++jd)
                        {
                            const int j_loc_shift = jd * db2;
                            const int j_bd_shift = j_mb_shift + j_loc_shift;
                            const int j_ub = min(j_bd_shift + db2, jb_ub) - j_bd_shift;
                            const int jd_column_w = min(db2, jb_column_w - j_loc_shift);
                            const int Abuf_diff_val = i_iter_count * jd_column_w;
                            double* Abuf_row = Abuf_dstripe +
                                id_block_h * j_loc_shift + i_lb * jd_column_w;

                            for (int kd = kd_lb; kd < kd_ub; ++kd)
                            {
                                const int k_loc_shift = kd * db1;
                                const int k_bd_shift = k_mb_shift + k_loc_shift;
                                const int k_lb = max(lambda, k_bd_shift) - k_bd_shift;
                                const int k_ub = min(k_bd_shift + db1, kb_ub) - k_bd_shift;
                                const int k_iter_count = k_ub - k_lb;
                                const int A_row_diff_val = k_iter_count * jd_column_w - 1;
                                const int kd_block_h = min(db1, kb_stripe_h - k_loc_shift);
                                const double* WY_row_shifted = WY_dstripe +
                                    id_block_h*k_loc_shift + i_lb*kd_block_h + k_lb;
                                const double* A_row = A_block +
                                    k_loc_shift * jb_column_w +
                                    kd_block_h * j_loc_shift + k_lb * jd_column_w;

                                for (int i = 0; i < i_iter_count; ++i)
                                {
                                    for (int j = 0; j < j_ub; ++j)
                                    {
                                        double sum(0.0);
                                        for (int k = 0; k < k_iter_count; ++k)
                                        {
                                            sum += WY_row_shifted[k] * (*A_row);
                                            A_row += jd_column_w;
                                        }
                                        Abuf_row[j] += sum;

                                        A_row -= A_row_diff_val;
                                    }

                                    WY_row_shifted += kd_block_h;
                                    Abuf_row += jd_column_w;
                                    A_row -= j_ub;
                                }
                                Abuf_row -= Abuf_diff_val;
                            }
                        }
                    }
                }
            }
        }

        copy_minor_double_block(A, Abuf, parameters, d1_shift, d2_shift);
    }// WHILE
    
    free(work_memory);

    return A;
}
