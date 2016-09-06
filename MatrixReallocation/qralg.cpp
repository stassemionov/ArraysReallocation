#include "qralg.h"
#include "service.h"

#include <algorithm>

using std::min;
using std::max;
using std::swap;

double* QR_WY_tiled(double* A, const TaskData& parameters)
{
    const int N = parameters.M_ROWS;
    const int b1 = parameters.B_ROWS;
    const int b2 = parameters.B_COLS;
    
    double* v = new double[N];
    double* w = new double[b2];
    double* W = new double[N * b2];
    double* Y = new double[N * b2];
    double* WY = new double[N * N];
    double* Abuf = new double[N * N];
    
    // lambda указывает на начало текущего блока
    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // Указывает на начало следующего блока за текущим
        const int t = min(lambda + b2, N);
        const int row_b_count = static_cast<int>(ceil(1.0 * (N - t) / b2));
        const int col_b_count = static_cast<int>(ceil(1.0 * (N - lambda) / b1));

        // Выполнение преобразований над текущим блоком
        for (int it = lambda; it < t; ++it)
        {
            memset(w, 0, (t - lambda)*sizeof(double));

            double norm = 0.0;
            double scalar = 1.0;
            double buf;
            // * Вычисление вектора Хаусхолдера
            for (int j = it; j < N; ++j)
            {
                buf = A[j * N + it];
                v[j] = buf;
                norm += buf * buf;
            }
            const double A_diag_el = A[it*N + it];
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
            for (int kb = 0; kb < col_b_count; ++kb)
            {
                const int lb = max(it, lambda + kb*b1);
                const int ub = min(lambda + (kb+1)*b1, N);
                for (int j = it; j < t; ++j)
                {
                    double sum = 0;
                    for (int k = lb; k < ub; ++k)
                    {
                        sum += A[k*N + j] * v[k];
                    }
                    w[j-lambda] += beta * sum;
                }
            }
            /*for (int k = it-lambda; k < b2; ++k)
            {
                w[k] *= beta;
            }*/

            // Преобразование текущего блока матрицы A
            // здесь блочность не нужна, она есть по умолчанию
            for (int j = it; j < N; ++j)
            {
                double* Aj = A + j * N;
                const double vj = v[j];
                for (int k = it; k < t; ++k)
                {
                    Aj[k] += vj * w[k - lambda];
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
        for (int it = 0; it < d; ++it)
        {
            const int shift = lambda + it;
            double scalar = 1, buf;
            // * Вычисление нормы и ск. произведения вектора Хаусхолдера
            for (int j = shift + 1; j < N; ++j)
            {
                buf = A[j * N + shift];
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
                    Y[i*b2] = v[i];
                    W[i*b2] = beta * v[i];
                }
            }
            else
            {
                memset(w, 0, d*sizeof(double));

                // Вычисление произведения (Y^t) * v
                for (int jb = 0; jb < col_b_count; ++jb)
                {
                    const int j_lb = max(shift, lambda+jb*b1);
                    const int j_ub = min(lambda+(jb+1)*b1, N);
                    for (int i = 0; i < it; ++i)
                    {
                        double sum = 0;
                        for (int j = j_lb; j < j_ub; ++j)
                        {
                            sum += Y[j*b2 + i] * v[j];
                        }
                        w[i] += sum;
                    }
                }

                // Вычисление произведения W * ((Y^t) * v)
                // здесь тоже блочность по умолчанию
                for (int i = lambda; i < N; ++i)
                {
                    double* Wi = W + i*b2;
                    double sum = 0;
                    for (int j = 0; j < it; ++j)
                    {
                        sum += Wi[j] * w[j];
                    }

                    buf = (i < shift) ? 0 : v[i];
                    W[i*b2 + it] = beta * (buf + sum);
                    Y[i*b2 + it] = buf;
                }
            }
        }

        // * Преобразование остальной части матрицы A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        // Умножение Y * W^t (матрицу WY используем повторно)
        for (int i = lambda; i < N; ++i)
        {
            const int kbound = min(d, i - lambda + 1);
            double* WYi = WY + i*N;
            const double* Yi = Y + i*b2;
            for (int j = lambda; j < N; ++j)
            {
                const double* Wj = W + j*b2;
                double sum = 0;
                for (int k = 0; k < kbound; ++k)
                {
                    sum += Yi[k] * Wj[k];
                }
                WYi[j] = sum;
            }
        }

        // Умножение A(lamda:M-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_lb = lambda + ib*b1;
            const int i_ub = min(i_lb + b1, N);
            
            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_lb = t + jb*b2;
                const int j_ub = min(j_lb + b2, N);

                // Для k = 0 отдельно, чтобы не выполнять проверку
                // в самом глубоком цикле
                const int k_lb0 = lambda;
                const int k_ub0 = min(k_lb0 + b1, N);
                for (int i = i_lb; i < i_ub; ++i)
                {
                    double* Abuf_i = Abuf + i*N;
                    double* Ai = A + i*N;
                    const double* WYi = WY + i*N;

                    for (int j = j_lb; j < j_ub; ++j)
                    {
                        double sum = Ai[j];
                        for (int k = k_lb0; k < k_ub0; ++k)
                        {
                            sum += WYi[k] * A[k*N + j];
                        }
                        Abuf_i[j] = sum;
                    }
                }

                for (int kb = 1; kb < col_b_count; ++kb)
                {
                    const int k_lb = lambda + kb*b1;
                    const int k_ub = min(k_lb + b1, N);

                    for (int i = i_lb; i < i_ub; ++i)
                    {
                        double* Abuf_i = Abuf + i*N;
                        const double* WYi = WY + i*N;

                        for (int j = j_lb; j < j_ub; ++j)
                        {
                            double sum = Abuf_i[j];
                            for (int k = k_lb; k < k_ub; ++k)
                            {
                                sum += WYi[k] * A[k*N + j];
                            }
                            Abuf_i[j] = sum;
                        }
                    }
                }
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

    delete[] v;
    delete[] w;
    delete[] W;
    delete[] Y;
    delete[] WY;
    delete[] Abuf;

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

    double* v = new double[N];
    double* w = new double[b2];
    double* W = new double[N * b2];
    double* Y = new double[N * b2];
    double* WY = new double[N * N];
    double* Abuf = new double[N * N];
    double* sum_vec = new double[db1];

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
        const int d1 = min(N, (curr_i_block_ind + 1)*b1) - d1_shift;
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
            double buf;
            const double* A_shifted_HV = A + full_shift * N + full_shift;
            double* v_shifted_HV = v + full_shift;
            for (int i = full_shift; i < N; ++i)
            {
                buf = *A_shifted_HV;
                (*v_shifted_HV++) = buf;
                norm += buf * buf;
                A_shifted_HV += N;
            }
            const double A_diag_el = A[full_shift*N + full_shift];
            const double A_diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + A_diag_el_sign * sqrt(norm);
            // используется нормировка, т.ч. v[it] = 1
            for (int j = full_shift + 1; j < N; ++j)
            {
                buf = (v[j] /= beta);
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
                    const int v_dif_val = i_ub - i_lb;
                    const double* A_dstripe = A + i_bd_shift * N;
                    const double* v_id_shifted = v + i_bd_shift + i_lb;

                    for (int jd = dblock_col_index; jd < dblock_count_in_curr_block_row; ++jd)
                    {
                        const int loc_drow_shift = jd * db2;
                        const int j_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int j_ub = min(d2 - loc_drow_shift, db2);
                        const int w_dif_val = j_ub - j_lb;
                        const int A_dif_val = N - w_dif_val;
                        double* w_shifted = w + loc_drow_shift + j_lb;
                        const double* A_row = (A_dstripe + i_lb*N) +
                            lambda + loc_drow_shift + j_lb;

                        for (int i = i_lb; i < i_ub; ++i)
                        {
                            const double vi = *v_id_shifted++;
                            for (int j = j_lb; j < j_ub; ++j)
                            {
                                (*w_shifted++) += (*A_row++) * vi;
                            }

                            w_shifted -= w_dif_val;
                            A_row += A_dif_val;
                        }
                        v_id_shifted -= v_dif_val;
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
                    const int v_diff_val = j_ub - j_lb;
                    const int j_lb_mult_N = j_lb * N;
                    double* A_dstripe = A + bd_shift * N;
                    const double* v_jd_shifted = v + bd_shift + j_lb;

                    for (int kd = dblock_col_index; kd < dblock_count_in_curr_block_row; ++kd)
                    {
                        const int loc_drow_shift = kd * db2;
                        const int k_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int dblock_w = min(db2, d2 - loc_drow_shift);
                        const int w_diff_val = dblock_w - k_lb;
                        const int A_diff_val = N - w_diff_val;
                        const double* w_shifted = w + loc_drow_shift + k_lb;
                        double* A_line = (A_dstripe + j_lb_mult_N) +
                            lambda + loc_drow_shift + k_lb;

                        for (int j = j_lb; j < j_ub; ++j)
                        {
                            const double vj = (*v_jd_shifted++);
                            for (int k = k_lb; k < dblock_w; ++k)
                            {
                                (*A_line++) += vj * (*w_shifted++);
                            }
                            A_line += A_diff_val;
                            w_shifted -= w_diff_val;
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
            double scalar = 1.0;
            double buf;
            double* A_shifted_HV = A + (full_shift + 1)*N + full_shift;
            double* v_shifted_HV = v + full_shift + 1;
            for (int j = full_shift + 1; j < N; ++j)
            {
                buf = *A_shifted_HV;
                scalar += buf * buf;
                (*v_shifted_HV++) = buf;
                A_shifted_HV += N;
            }
            v[full_shift] = 1.0;
            double beta = -2.0 / scalar;

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
                double buf_WYinit;
                for (int i = lambda; i < N; ++i)
                {
                    buf_WYinit = (*v_shifted_WYinit++);
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

                            for (int i = i_lb; i < i_ub; ++i)
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
                        const int i_lb_mult_b2 = i_lb * b2;
                        double* Wdstripe = W + bd_shift*b2;
                        double* Ydstripe = Y + (Wdstripe - W);
                        // Vector to contain results of production of W's row and vector w
                        memset(sum_vec, 0, sizeof(double)*db1);

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it, loc_drow_shift + db2) - loc_drow_shift;
                            const int diff_val = b2 - j_ub;
                            const double* w_shifted = w + loc_drow_shift;
                            const double* Wrow = Wdstripe + i_lb_mult_b2 + loc_drow_shift;
                            double* sum_vec_shifted = sum_vec + i_lb;

                            for (int i = i_lb; i < i_ub; ++i)
                            {
                                double sum_loc = 0.0;
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
                        for (int i = i_lb; i < i_ub; ++i)
                        {
                            buf = (*v_shifted_WYadd++);
                            *Wdblock_add = beta * (buf + (*sum_vec_shifted++));
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
        for (int i = d1_shift; i < N; ++i)
        {
            memset(WY_shifted, 0, (N - d1_shift)*sizeof(double));
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
                            for (int i = i_lb; i < i_ub; ++i)
                            {
                                const double* const Yrow_save = Yrow;
                                const double* Wrow = (Wdstripe + j_lb_mul_b2) + k_loc_shift;
                                for (int j = j_lb; j < j_ub; ++j)
                                {
                                    double sum = 0.0;
                                    for (int k = 0; k < k_ub; ++k)
                                    {
                                        sum += (*Yrow++) * (*Wrow++);
                                    }
                                    (*WYrow++) += sum;

                                    Yrow = Yrow_save;
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
        for (int i = d1_shift; i < N; ++i)
        {
            memcpy(
                Abuf_shifted,
                A_shifted,
                (N - d2_shift) * sizeof(double));
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
                                // Points to 'i_lb'-th row of ['id','kd'] small tile of 'WY'
                                const double* WY_dstripe_shifted = (WY_dstripe + i_lb_mul_N) + k_bd_shift;
                                // Points to first row of ['kd','jd'] small tile of 'A'
                                const double* A_dstripe_shifted = (A + k_bd_shift*N) + j_bd_shift;

                                double* Abuf_row = Abuf_dstripe_shifted;
                                for (int i = i_lb; i < i_ub; ++i)
                                {
                                    const double* WY_row = WY_dstripe_shifted;
                                    for (int j = 0; j < j_ub; ++j)
                                    {
                                        const double* A_row = A_dstripe_shifted + j;
                                        double sum = 0.0;
                                        for (int k = 0; k < k_ub; ++k)
                                        {
                                            sum += (*WY_row++) * (*A_row);
                                            A_row += N;
                                        }
                                        (*Abuf_row++) += sum;
                                        
                                        WY_row -= k_ub;
                                    }
                                    Abuf_row += Abuf_dif_value;
                                    WY_dstripe_shifted += N;
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
            memcpy(
                A_shifted_toUpdate,
                Abuf_shifted_toUpdate,
                (N - d2_shift) * sizeof(double));
            A_shifted_toUpdate += N;
            Abuf_shifted_toUpdate += N;
        }
    }// WHILE

    delete[] v;
    delete[] w;
    delete[] W;
    delete[] Y;
    delete[] WY;
    delete[] Abuf;
    delete[] sum_vec;

    return A;
}

double* QR_WY_standard(double* A, const int N, const int r)
{
    double* v    = new double[N];
    double* w    = new double[r];
    double* W    = new double[N * r];
    double* Y    = new double[N * r];
    double* WY   = new double[N * N];
    double* Abuf = new double[N * N];

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

    delete[] v;
    delete[] w;
    delete[] W;
    delete[] Y;
    delete[] WY;
    delete[] Abuf;

    return A;
}

double* QR_WY_block(double* A, const TaskData& parameters)
{
    const int N = parameters.M_ROWS;
    const int b1 = parameters.B_ROWS;
    const int b2 = parameters.B_COLS;
    const int block_count_in_col = static_cast<int>(ceil(1.0 * N / b1));
    const int block_count_in_row = static_cast<int>(ceil(1.0 * N / b2));

    double* v = new double[N];
    double* w = new double[b2];
    double* W = new double[N*b2];
    double* Y = new double[N*b2];
    double* WY = new double[N*N];
    double* Abuf = new double[N*N];

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
        // сдвиги по вертикали и горизонтали для текущего блока
        const int d1_shift = curr_i_block_ind*b1;
        const int d2_shift = curr_j_block_ind*b2;

        // действительные величины высоты/ширины текущих блочных полосы/столбца
        const int d1 = min(N, (curr_i_block_ind + 1)*b1) - d1_shift;
        const int d2 = t - lambda;
        
        // * Выполнение преобразований над текущим блочным столбцом
        for (int it = 0; it < d2; ++it)
        {
            // сдвиг до 'lambda+it'-го элемнта в столбце
            const int full_shift = lambda + it;
            // фиксация перехода в следующий по вертикали блок
            const int b_begin = (full_shift - d1_shift) / b1;

            memset(w, 0, d2*sizeof(double));
            memset(v + lambda, 0, it*sizeof(double));

            double norm = 0;
            double scalar = 1;
            double buf;

            // * Вычисление вектора Хаусхолдера
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                // сдвиг линии блока от линии матрицы
                const int mb_shift = jb*b1 + d1_shift;
                // верхняя и нижняя границы чтения строк в блоке jb (локальная нумерация)
                const int j_lb = max(full_shift, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift+b1, N) - mb_shift;
                // высота блока jb
                const int block_h = min(b1, N - mb_shift);
                // указатель на начало блока jb
                const double* Ablock = A + mb_shift*N + block_h*d2_shift;
                double* v_jb_shift = v + mb_shift;
                // j - номер строки в блоке jb (локальная нумерация)
                for (int j = j_lb; j < j_ub; ++j)
                {
                    buf = Ablock[j*d2 + it];
                    v_jb_shift[j] = buf;
                    norm += buf * buf;
                }
            }

            const int diag_el_block_h = min( b1, N - (d1_shift + b_begin*b1) );
            const double A_diag_el = A[d1_shift*N + b_begin*b1*N + d2_shift*diag_el_block_h + (full_shift % b1)*d2 + it];
            const double diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + diag_el_sign * sqrt(norm);
            // используется нормировка, т.ч. v[it] = 1
            for (int j = full_shift + 1; j < N; ++j)
            {
                buf = (v[j] /= beta);
                scalar += buf * buf;
            }
            v[full_shift] = 1.0;
            beta = -2.0 / scalar;
        
            // * Применение преобразования
            // Вычисление вектора w = beta * (A(:,lambda:t-1)^t) * v .
            // начинаем с it, т.к. it-ая матрица Хаусхолдера действует на матрицу A(it:,it:)
            for (int kb = b_begin; kb < col_b_count; ++kb)
            {
                // сдвиг линии блока от линии матрицы
                const int mb_shift = kb*b1 + d1_shift;
                // верхняя и нижняя границы чтения строк в блоке jb (локальная нумерация)
                const int k_lb = max(full_shift, mb_shift) - mb_shift;
                const int k_ub = min(mb_shift + b1, N) - mb_shift;
                // высота блока jb
                const int block_h = min(b1, N - mb_shift);
                // указатель на начало блока jb
                const double* Ablock = A + mb_shift*N + block_h*d2_shift;
                const double* v_kb_shift = v + mb_shift;
                // j - номер столбца в блоке kb (локальная нумерация)
                for (int j = it; j < d2; ++j)
                {
                    double sum = 0;
                    for (int k = k_lb; k < k_ub; ++k)
                    {
                        sum += Ablock[k*d2 + j] * v_kb_shift[k];
                    }
                    w[j] += sum;
                }
            }
            for (int k = it; k < d2; ++k)
            {
                w[k] *= beta;
            }

            // Преобразование текущего блока матрицы A
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                // сдвиг линии блока от линии матрицы
                const int mb_shift = d1_shift + jb*b1;
                // верхняя и нижняя границы чтения строк в блоке jb (локальная нумерация)
                const int j_lb = max(full_shift, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // высота блока jb
                const int j_block_h = min(b1, N - mb_shift);
                // указатель на начало блока jb
                double* Ablock = A + mb_shift*N + j_block_h*d2_shift;
                // j - номер строки в блоке jb (локальная нумерация)
                for (int j = j_lb; j < j_ub; ++j)
                {
                    // указатель на начало j-ой строки блока jb
                    double* Aj = Ablock + j*d2;
                    const double vj = v[mb_shift + j];
                    // k - номер столбца в блоке jb (локальная нумерация)
                    for (int k = it; k < d2; ++k)
                    {
                        Aj[k] += vj * w[k];
                    }
                }
            }
            
            // * Заполнение поддиагнальной части it-го столбца матрицы A
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
                double* Ablock = A + mb_shift*N + j_block_h*d2_shift;
                const double* v_shift = v + mb_shift;
                // j - номер строки в блоке jb (локальная нумерация)
                for (int j = j_lb; j < j_ub; ++j)
                {
                    Ablock[j*d2 + it] = v_shift[j];
                }
            }
        }// FOR


        memset(W, 0, N*b2*sizeof(double));
        memset(Y, 0, N*b2*sizeof(double));

        // * Вычисление WY-представления произведения матриц Хаусхолдера
        for (int it = 0; it < d2; ++it)
        {
            // сдвиг до 'lambda+it'-го элемнта в столбце
            const int full_shift = lambda + it;
            // фиксация перехода в следующий по вертикали блок
            const int b_begin = (full_shift - d1_shift) / b1;

            memset(w, 0, d2*sizeof(double));
            memset(v + lambda, 0, it*sizeof(double));

            double scalar = 1, buf;
            // * Вычисление нормы и ск. произведения вектора Хаусхолдера
            for (int jb = (full_shift + 1 - d1_shift) / b1; jb < col_b_count; ++jb)
            {
                // сдвиг линии блока от линии матрицы
                const int mb_shift = jb*b1 + d1_shift;
                // верхняя и нижняя границы чтения строк в блоке jb (локальная нумерация)
                const int j_lb = max(full_shift+1, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // высота блока jb
                const int block_h = min(b1, N - mb_shift);
                // указатель на начало блока jb
                const double* Ablock = A + mb_shift*N + block_h*d2_shift;
                double* v_jb_shift = v + mb_shift;
                // j - номер строки в блоке jb (локальная нумерация)
                for (int j = j_lb; j < j_ub; ++j)
                {
                    buf = Ablock[j*d2 + it];
                    scalar += buf * buf;
                    v_jb_shift[j] = buf;
                }
            }
            v[full_shift] = 1.0;
            double beta = -2.0 / scalar;
            
            // it - число накопленных столбцов в матрицах W и Y (= номеру дозаписываемого столбца)
            if (it == 0)
            {
                // Начальное заполнение первого столбца
                // замечание: первые lambda строк матриц W и Y - нулевые,
                //            матрицы имеют ступенчатый вид.
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    // сдвиг линии блока от линии матрицы
                    const int mb_shift = ib*b1 + d1_shift;
                    // верхняя и нижняя границы чтения строк в блоке ib (локальная нумерация)
                    const int i_lb = max(full_shift, mb_shift) - mb_shift;
                    const int i_ub = min(mb_shift + b1, N) - mb_shift;
                    // высота блока ib
                    const int block_h = min(b1, N - mb_shift);
                    // указатели на начала блоков ib матриц Y и W
                    double* Yblock = Y + mb_shift*b2;
                    double* Wblock = W + mb_shift*b2;
                    const double* v_shift = v + mb_shift;
                    for (int i = i_lb; i < i_ub; ++i)
                    {
                        Yblock[i*b2] = v_shift[i];
                        Wblock[i*b2] = beta * v_shift[i];
                    }
                }
            }
            else
            {
                // Вычисление произведения (Y^t) * v
                for (int jb = b_begin; jb < col_b_count; ++jb)
                {
                    // сдвиг линии блока от линии матрицы
                    const int mb_shift = jb*b1 + d1_shift;
                    // верхняя и нижняя границы чтения строк в блоке jb (локальная нумерация)
                    const int j_lb = max(full_shift, mb_shift) - mb_shift;
                    const int j_ub = min(mb_shift + b1, N) - mb_shift;
                    // высота блока jb
                    const int block_h = min(b1, N - mb_shift);
                    // указатель на начало блока jb
                    const double* Yblock = Y + mb_shift*b2;
                    double* v_shift = v + mb_shift;
                    // j - номер строки в блоке jb (локальная нумерация)
                    for (int j = j_lb; j < j_ub; ++j)
                    {
                        const double* Yline = Yblock + j*b2;
                        const double vj = v_shift[j];
                        for (int i = 0; i < it; ++i)
                        {
                            w[i] += Yline[i] * vj;
                        }
                    }
                }
                
                // Вычисление произведения W * ((Y^t) * v)
                // здесь тоже блочность по умолчанию
                for (int ib = 0; ib < col_b_count; ++ib)
                {
                    // сдвиг линии блока от линии матрицы
                    const int mb_shift = ib*b1 + d1_shift;
                    // верхняя и нижняя границы чтения строк в блоке ib (локальная нумерация)
                    const int i_lb = max(lambda, mb_shift) - mb_shift;
                    const int i_ub = min(b1, N - mb_shift);
                    // высота блока ib
                    const int block_h = min(b1, N - mb_shift);
                    // указатели на начала блоков ib матриц Y и W
                    double* Yblock = Y + mb_shift*b2;
                    double* Wblock = W + (Yblock - Y);
                    const double* v_shift = v + mb_shift;
                    for (int i = i_lb; i < i_ub; ++i)
                    {
                        double* Wline = Wblock + i*b2;
                        double sum = 0;
                        for (int j = 0; j < it; ++j)
                        {
                            sum += Wline[j] * w[j];
                        }
                        buf = v_shift[i];
                        Wline[it] = beta * (buf + sum);
                        Yblock[i*b2+it] = buf;
                    }
                }
            }
        }

        // * Преобразование остальной части матрицы A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)
        
        // Умножение Y * W^t (матрицу WY используем повторно).
        // Матрица WY будет иметь b1 x b1 размещение
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int i_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int i_ub = min(b1, N - i_mb_shift);
            const double* Yblock = Y + i_mb_shift*b2;

            const int ib_stripe_h = min(b1, N - i_mb_shift);
            double* WYstripe = WY + i_mb_shift*N + ib_stripe_h*d1_shift;

            for (int jb = 0; jb < col_b_count; ++jb)
            {
                const int j_mb_shift = d1_shift + jb*b1;
                const int j_lb = max(lambda, j_mb_shift) - j_mb_shift;
                const int j_ub = min(b1, N - j_mb_shift);
                const double* Wblock = W + j_mb_shift*b2;

                const int jb_block_w = min(b1, N - j_mb_shift);
                double* WYblock = WYstripe + jb*ib_stripe_h*b1;

                for (int i = i_lb; i < i_ub; ++i)
                {
                    // ограничение, связанное со ступенчатым видом матрицы Y
                    const int kbound = min(t, i_mb_shift + i + 1) - lambda;

                    const double* Ybi = Yblock + i*b2;
                    double* WYb_line = WYblock + i*jb_block_w;
                                        
                    for (int j = j_lb; j < j_ub; ++j)
                    {
                        const double* Wbj = Wblock + j*b2;
                        double sum = 0.0;
                        for (int k = 0; k < kbound; ++k)
                        {
                            sum += Ybi[k] * Wbj[k];
                        }
                        WYb_line[j] = sum;
                    }
                }
            }
        }

        // Умножение A(lamda:M-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1)
        // Размещение: (b1 x b2) = (b1 x b1) * (b1 x b2)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int i_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int i_ub = min(b1, N - i_mb_shift);
            const int ib_stripe_h = min(b1, N - i_mb_shift);

            // указатели сдвинуты до левой границы записываемой/читаемой области
            const double* Agen_ib_stripe_shifted = A + i_mb_shift*N + t*ib_stripe_h;
            double* Abuf_ib_stripe_shifted = Abuf + (Agen_ib_stripe_shifted - A);
            const double* WY_ib_stripe_shifted = WY + i_mb_shift*N + d1_shift*ib_stripe_h;

            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_mb_shift = t + jb*b2;
                const int jb_block_w = min(b2, N - j_mb_shift);

                // этот блок нужен для чтения на первой итерации по kb
                const double* Agen_init_ib_jb_block = Agen_ib_stripe_shifted + jb*ib_stripe_h*b2;
                // в этот блок будем записывать результат произведения
                double* Abuf_ib_jb_block = Abuf_ib_stripe_shifted +
                    (Agen_init_ib_jb_block - Agen_ib_stripe_shifted);

                // Для kb = 0 отдельно, чтобы не выполнять проверку
                // для присваивания в самом глубоком цикле (в А или в Abuf)
                const int k0_lb = max(lambda, d1_shift) - d1_shift;
                const int kb0_block_size = min(b1, N - d1_shift);

                const double* Aright_block_k0_jb = A + d1_shift*N + j_mb_shift*kb0_block_size;

                for (int i = i_lb; i < i_ub; ++i)
                {
                    double* Abuf_i = Abuf_ib_jb_block + i*jb_block_w;
                    const double* Agen_init_i = Agen_init_ib_jb_block + (Abuf_i - Abuf_ib_jb_block);
                    const double* WYi = WY_ib_stripe_shifted + i*kb0_block_size;

                    for (int j = 0; j < jb_block_w; ++j)
                    {
                        double sum = Agen_init_i[j];
                        const double* Aright_block_k0_jb_shifted = Aright_block_k0_jb +
                            k0_lb*jb_block_w + j;
                        for (int k = k0_lb; k < kb0_block_size; ++k)
                        {
                            sum += WYi[k] * (*Aright_block_k0_jb_shifted);
                            Aright_block_k0_jb_shifted += jb_block_w;
                        }
                        Abuf_i[j] = sum;
                    }
                }

                for (int kb = 1; kb < col_b_count; ++kb)
                {
                    const int k_mb_shift = d1_shift + kb*b1;
                    const int kb_block_size = min(b1, N - k_mb_shift);
                    
                    const double* WYblock = WY_ib_stripe_shifted + kb*ib_stripe_h*b1;
                    const double* Ablock = A + k_mb_shift*N + j_mb_shift*kb_block_size;
                    for (int i = i_lb; i < i_ub; ++i)
                    {
                        double* Abuf_i = Abuf_ib_jb_block + i*jb_block_w;
                        const double* WYi = WYblock + i*kb_block_size;
                        for (int j = 0; j < jb_block_w; ++j)
                        {
                            double sum = 0.0;
                            const double* Ablock_shifted = Ablock + j;
                            for (int k = 0; k < kb_block_size; ++k)
                            {
                                sum += WYi[k] * (*Ablock_shifted);
                                Ablock_shifted += jb_block_w;
                            }
                            Abuf_i[j] += sum;
                        }
                    }
                }
            }
        }
        
        // * копирование преобразованной части матрицы A
        // копирование первой блочное полосы (отдельно, т.к. она может быть непоной по высоте)
        const int i_lb = max(lambda, d1_shift) - d1_shift;
        const int i_ub = min(d1_shift + b1, N) - d1_shift;
        const int ib0_stripe_h = min(b1, N - d1_shift);

        const double* Abuf_stripe = Abuf + d1_shift*N;
        double* Astripe = A + d1_shift*N;
        for (int jb = 0; jb < row_b_count; ++jb)
        {
            const int j_mb_shift = t + jb*b2;
            const int jb_block_w = min(b2, N - j_mb_shift);

            double* A_block_shifted = Astripe + ib0_stripe_h*j_mb_shift + i_lb*jb_block_w;
            const double* Abuf_block_shifted = Abuf_stripe + ib0_stripe_h*j_mb_shift + i_lb*jb_block_w;
         
            memcpy(A_block_shifted, Abuf_block_shifted, (i_ub - i_lb)*jb_block_w*sizeof(double));
        }

        for (int ib = 1; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int ib_stripe_h = min(b1, N - i_mb_shift);

            const double* Abuf_stripe_shifted = Abuf + i_mb_shift*N + t*ib_stripe_h;
            double* Astripe_shifted = A + i_mb_shift*N + t*ib_stripe_h;

            memcpy(Astripe_shifted, Abuf_stripe_shifted, ib_stripe_h*(N - t)*sizeof(double));
        }
    }// WHILE

    delete[] v;
    delete[] w;
    delete[] W;
    delete[] Y;
    delete[] WY;
    delete[] Abuf;

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

    double* v = new double[N];
    double* w = new double[b2];
    double* W = new double[N * b2];
    double* Y = new double[N * b2];
    double* WY = new double[N * N];
    double* Abuf = new double[N * N];
    double* sum_vec = new double[db1];

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
        // Real values of height/width of current block-stripe/column
        const int d1 = min(N, (curr_i_block_ind + 1)*b1) - d1_shift;
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
            
            double norm = 0.0;
            double scalar = 1.0;
            double buf;

            // * Householder vector computation...
            // Using of 'it'-th column of current block-column,
            // starting with (lambda+it)-th position
            for (int ib = b_begin; ib < col_b_count; ++ib)
            {
                // Shift of current big block first line
                // from matrix first line
                const int mb_shift = d1_shift + ib*b1;
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
                const double* Ablock = A + mb_shift*N + ib_block_h*d2_shift;

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
                    const int bd_shift = loc_shift + mb_shift;
                    // Bounds for current small block reading.
                    // Reasons of using are the same as 'i_lb' and 'i_ub'
                    const int lb = max(full_shift, bd_shift) - bd_shift;
                    const int ub = min(bd_shift + db1, i_ub) - bd_shift;
                    // 'id'-th small block height (number of rows)
                    const int id_block_h = min(db1, ib_block_h - loc_shift);
                    // Pointer to beginning of required small block in 'id'-th stripe
                    const double* Adblock_shifted = Ablock + loc_shift*d2 +
                        dblock_col_index*id_block_h*db2 + loc_it;
                    double* v_shifted = v + bd_shift + lb;
                    
                    for (int i = lb*dblock_width; i < ub*dblock_width; i += dblock_width)
                    {
                        buf = Adblock_shifted[i];
                        *(v_shifted++) = buf;
                        norm += buf * buf;
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
            double beta = 1.0 / (A_diag_el + diag_el_sign * sqrt(norm));
            // Use normalization such that v[it] = 1
            for (int j = full_shift + 1; j < N; ++j)
            {
                buf = (v[j] *= beta);
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
                        // Pointer to beginning of 'jb'-th small block
                        const double* Adblock = Adstripe + dblock_h*loc_drow_shift;
                        double* w_shifted = w + loc_drow_shift;

                        const int klb = lb * dblock_w;
                        const int kub = ub * dblock_w;
                        for (int j = j_lb; j < j_ub; ++j)
                        {
                            double sum = 0.0;
                            const double* Adblock_shifted = Adblock + j;
                            for (int k = klb; k < kub; k += dblock_w)
                            {
                                sum += Adblock_shifted[k] * (*(v_kd_shifted++));
                            }
                            v_kd_shifted -= (ub - lb);
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
                const int mb_shift = d1_shift + jb*b1;
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
                    const int jd_block_h = min(db1, jb_block_h - loc_shift);
                    double* Adstripe = Ablock + loc_shift*d2;
                    const double* v_shifted = v + bd_shift;

                    for (int kd = dblock_col_index; kd < dblock_count_in_curr_block_row; ++kd)
                    {
                        const int loc_drow_shift = kd * db2;
                        const int k_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int dblock_w = min(db2, d2 - loc_drow_shift);
                        double* Adblock_shifted = Adstripe + jd_block_h*loc_drow_shift +
                            k_lb + j_lb*dblock_w;
                        const double* w_shifted = w + loc_drow_shift;

                        for (int j = j_lb; j < j_ub; ++j)
                        {
                            double* Aj = Adblock_shifted;
                            Adblock_shifted += dblock_w;
                            const double vj = v_shifted[j];
                            for (int k = k_lb; k < dblock_w; ++k)
                            {
                                (*Aj++) += vj * w_shifted[k];
                            }
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
                    double* Adblock_shifted = Ablock + loc_shift*d2 +
                        dblock_col_index*id_block_h*db2 + loc_it;
                    const double* v_shifted = v + bd_shift + i_lb;
                    
                    for (int i = i_lb*dblock_width; i < i_ub*dblock_width; i += dblock_width)
                    {
                        Adblock_shifted[i] = *(v_shifted++);
                    }
                }
            }
        }// FOR

        // * Computation of WY-representation of Householder matrices production...
        memset(W + d1_shift*b2, 0, min(bsizes_ratio*b1, N - d1_shift)*b2*sizeof(double));
        memset(Y + d1_shift*b2, 0, min(bsizes_ratio*b1, N - d1_shift)*b2*sizeof(double));
        for (int it = 0; it < d2; ++it)
        {
            memset(w, 0, d2*sizeof(double));
            memset(v + lambda, 0, it*sizeof(double));

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
            double buf;
            for (int ib = (full_shift+1 - d1_shift)/b1; ib < col_b_count; ++ib)
            {
                const int mb_shift = d1_shift + ib*b1;
                const int i_lb = max(full_shift + 1, mb_shift) - mb_shift;
                const int i_ub = min(mb_shift + b1, N);
                const int ib_block_h = min(b1, N - mb_shift);
                const double* Ablock = A + mb_shift*N + ib_block_h*d2_shift;

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

                    for (int i = lb; i < ub; ++i)
                    {
                        buf = *Adblock_shifted;
                        Adblock_shifted += dblock_width;
                        scalar += buf * buf;
                        *(v_shifted++) = buf;
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
                    const int ib_block_h = min(b1, N - mb_shift);

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
                        const int id_block_h = min(db1, ib_block_h - loc_shift);
                        // Width of block is 'b2',
                        // because memory was allocated for complete matrices
                        double* Ydblock = Y + bd_shift*b2 + lb*dblock_width_WY;
                        double* Wdblock = W + (Ydblock - Y);
                        double* v_shifted = v + bd_shift + lb;

                        for (int i = lb; i < ub; ++i)
                        {
                            *Ydblock = *v_shifted;
                            *Wdblock = beta * (*v_shifted++);
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
                        const int id_block_h = min(db1, ib_block_h - loc_shift);
                        const double* Ydstripe = Y + bd_shift*b2;
                        const double* v_shifted = v + bd_shift;

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it - loc_drow_shift, db2);
                            const int dblock_w = min(db2, b2 - loc_drow_shift);
                            const double* Ydblock = Ydstripe + id_block_h*loc_drow_shift;
                            double* w_shifted = w + loc_drow_shift;

                            for (int i = i_lb; i < i_ub; ++i)
                            {
                                const double* Yline = Ydblock + i*dblock_w;
                                const double vi = v_shifted[i];
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    w_shifted[j] += Yline[j] * vi;
                                }
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
                        const int id_block_h = min(db1, ib_block_h - loc_shift);
                        double* Wdstripe = W + bd_shift*b2;
                        double* Ydstripe = Y + (Wdstripe - W);
                        // Vector to contain results of production of W's row and vector w
                        memset(sum_vec, 0, db1 << 3);

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it, loc_drow_shift + db2) - loc_drow_shift;
                            const int dblock_w = min(db2, b2 - loc_drow_shift);
                            double* Wdblock = Wdstripe + id_block_h*loc_drow_shift;
                            const double* w_shifted = w + loc_drow_shift;

                            for (int i = i_lb; i < i_ub; ++i)
                            {
                                double* Wline = Wdblock + i*dblock_w;
                                double sum_loc = 0.0;
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    sum_loc += Wline[j] * w_shifted[j];
                                }
                                // Scalar production of W matrix line and vector w
                                sum_vec[i] += sum_loc;
                            }
                        }
                        // Insertion of computed values into 'it'-th column of Y and W
                        const double* v_shifted = v + bd_shift + i_lb;
                        // Pointers to small blocks containing 'it'-th column of Y and W
                        double* Wdblock_add = Wdstripe + dblock_col_index*id_block_h*db2 + loc_it;
                        double* Ydblock_add = Ydstripe + (Wdblock_add - Wdstripe);
                        const double* sum_vec_shifted = sum_vec + i_lb;
                        for (int i = i_lb*dblock_width_WY; i < i_ub*dblock_width_WY; i += dblock_width_WY)
                        {
                            buf = (*v_shifted++);
                            Wdblock_add[i] = beta * (buf + (*sum_vec_shifted++));
                            Ydblock_add[i] = buf;
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
        filmat_part(WY, 0, WY_task_class.getDataRef(), d1_shift, d1_shift);

        // Multiplication Y * W^t .
        // Remark: WY matrix get (B1,B1,D1,D1)-layout.
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
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
                const int j_mb_shift = d1_shift + jb*b1;
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
                        const int j_ub = min(db1, jb_ub - j_bd_shift) ;
                        const int jd_block_h = min(db1, jb_stripe_h - j_loc_shift);
                        const double* Wdstripe = W + j_bd_shift*b2;
                        double* WYdblock = WYdstripe + id_block_h*j_loc_shift;
                        
                        for (int kd = 0; kd < kd_ud; ++kd)
                        {
                            const int k_loc_shift = kd * db2;
                            const int k_ub = min(d2 - k_loc_shift, db2);
                            const int kd_dblock_w = min(db2, b2 - k_loc_shift);
                            const double* Ydblock = Ydstripe + id_block_h*k_loc_shift;
                            const double* Wdblock = Wdstripe + jd_block_h*k_loc_shift;

                            for (int i = i_lb; i < i_ub; ++i)
                            {
                                double* WYrow = WYdblock + i*jd_block_h + j_lb;
                                const double* Yrow = Ydblock + i*kd_dblock_w;
                                
                                for (int j = j_lb; j < j_ub; ++j)
                                {
                                    const double* Wrow = Wdblock + j*kd_dblock_w;
                                    double sum = 0.0;
                                    for (int k = 0; k < k_ub; ++k)
                                    {
                                        sum += Yrow[k] * (*Wrow++);
                                    }
                                    (*WYrow++) += sum;
                                }
                            }
                        }
                    }
                }
            }
        }

        filmat_part(Abuf, A, parameters, d1_shift, d2_shift);

        // Multiplication A(lamda:N-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1) <=>
        //                A(lamda:N-1, t:N-1) += WY * A(lamda:N-1, t:N-1) .
        // Layout pattern: (b1 x b2, d1 x d2) = (b1 x b1, d1 x d1) * (b1 x b2, d1 x d2).
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int ib_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int ib_ub = min(i_mb_shift + b1, N);
            const int ib_stripe_h = min(b1, N - i_mb_shift);
            // Pointer to beginning of 'ib'-th big-block-stripe of WY matrix
            const double* WY_stripe = WY + i_mb_shift*N;
            // Pointer to beginning of 'ib'-th big-block-stripe of buffer matrix
            double* Abuf_stripe = Abuf + (WY_stripe - WY);

            const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
            const int id_ub = (ib == col_b_count - 1) ?
                dblock_count_in_diff_bcol :
                dblock_count_in_bcol;
            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_mb_shift = t + jb*b2;
                const int jb_ub = min(j_mb_shift + b2, N);
                const int jb_column_w = min(b2, N - j_mb_shift);
                double* Abuf_block = Abuf_stripe + j_mb_shift*ib_stripe_h;

                const int jd_ub = (jb == row_b_count - 1) ?
                    dblock_count_in_diff_brow :
                    dblock_count_in_brow;
                for (int kb = 0; kb < col_b_count; ++kb)
                {
                    const int k_mb_shift = d1_shift + kb*b1;
                    const int kb_lb = max(lambda, k_mb_shift) - k_mb_shift;
                    const int kb_ub = min(k_mb_shift + b1, N);
                    const int kb_stripe_h = min(b1, N - k_mb_shift);
                    const double* WY_block = WY_stripe + ib_stripe_h*k_mb_shift;
                    const double* A_block = A + k_mb_shift*N + j_mb_shift*kb_stripe_h;

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
                        const int id_block_h = min(db1, ib_stripe_h - i_loc_shift);
                        const double* WY_dstripe = WY_block + i_loc_shift*kb_stripe_h;
                        double* Abuf_dstripe = Abuf_block + i_loc_shift*jb_column_w;

                        for (int jd = 0; jd < jd_ub; ++jd)
                        {
                            const int j_loc_shift = jd * db2;
                            const int j_bd_shift = j_mb_shift + j_loc_shift;
                            const int j_ub = min(j_bd_shift + db2, jb_ub) - j_bd_shift;
                            const int jd_column_w = min(db2, jb_column_w - j_loc_shift);
                            double* Abuf_dblock = Abuf_dstripe + id_block_h*j_loc_shift;

                            // Remaining 'kd'-loop iterations
                            for (int kd = kd_lb; kd < kd_ub; ++kd)
                            {
                                const int k_loc_shift = kd * db1;
                                const int k_bd_shift = k_mb_shift + k_loc_shift;
                                const int k_ub = min(k_bd_shift + db1, kb_ub) - k_bd_shift;
                                const int kd_block_h = min(db1, kb_stripe_h - k_loc_shift);
                                const double* WY_dblock = WY_dstripe + id_block_h*k_loc_shift;
                                const double* A_dblock = A_block +
                                    k_loc_shift*jb_column_w + kd_block_h*j_loc_shift;

                                for (int i = i_lb; i < i_ub; ++i)
                                {
                                    double* Abuf_row = Abuf_dblock + i*jd_column_w;
                                    const double* WY_row = WY_dblock + i*kd_block_h;
                                    for (int j = 0; j < j_ub; ++j)
                                    {
                                        const double* A_dblock_shifted = A_dblock + j;
                                        double sum = 0.0;
                                        for (int k = 0; k < k_ub; ++k)
                                        {
                                            sum += WY_row[k] * (*A_dblock_shifted);
                                            A_dblock_shifted += jd_column_w;
                                        }
                                        (*Abuf_row++) += sum;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        filmat_part(A, Abuf, parameters, d1_shift, d2_shift);
    }// WHILE
    
    delete[] v;
    delete[] w;
    delete[] W;
    delete[] Y;
    delete[] WY;
    delete[] Abuf;
    delete[] sum_vec;

    return A;
}
