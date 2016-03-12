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
    double* w = new double[N];
    double* W = new double[N * b2];
    double* Y = new double[N * b2];
    double* WY = new double[N * N];
    double* Abuf = new double[N * N];
    
    // lambda указывает на начало текущего блока
    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // Указывает на начало следующего блока за текущим
        const int t = min(lambda + b2, N);
        const int alloc_size = (N - lambda)*sizeof(double);
        const int row_b_count = static_cast<int>(ceil(1.0 * (N - t) / b2));
        const int col_b_count = static_cast<int>(ceil(1.0 * (N - lambda) / b1));

        // Выполнение преобразований над текущим блоком
        for (int it = lambda; it < t; ++it)
        {
            // ничего, если занулятся несколько лишних коодинат
            memset(w + lambda, 0, alloc_size);

            double norm = 0;
            double scalar = 1;
            double buf;
            // * Вычисление вектора Хаусхолдера
            for (int j = it; j < N; ++j)
            {
                buf = A[j * N + it];
                v[j] = buf;
                norm += buf * buf;
            }
            const double A_diag_el = A[it * N + it];
            double beta = A_diag_el + A_diag_el / abs(A_diag_el) * sqrt(norm);
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
                    w[j] += beta * sum;
                }
            }
     
            // Преобразование текущего блока матрицы A
            // здесь блочность не нужна, она есть по умолчанию
            for (int j = it; j < N; ++j)
            {
                double* Aj = A + j * N;
                const double vj = v[j];
                for (int k = it; k < t; ++k)
                {
                    Aj[k] += vj * w[k];
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
                const int k_lb = lambda;
                const int k_ub = min(k_lb + b1, N);
                for (int i = i_lb; i < i_ub; ++i)
                {
                    double* Abuf_i = Abuf + i*N;
                    double* Ai = A + i*N;
                    const double* WYi = WY + i*N;

                    for (int j = j_lb; j < j_ub; ++j)
                    {
                        double sum = Ai[j];
                        for (int k = k_lb; k < k_ub; ++k)
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
                        double* Ai = A + i*N;
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

double* QR_WY_standard(double* A, const int N, const int r)
{
    double* v    = new double[N];
    double* w    = new double[N];
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
            double norm = 0;
            double scalar = 1;
            double buf;
            // * Вычисление вектора Хаусхолдера
            for (int j = it; j < N; ++j)
            {
                buf = A[j * N + it];
                v[j] = buf;
                norm += buf * buf;
            }
            const double A_diag_el = A[it * N + it];
            double beta = A_diag_el + A_diag_el / abs(A_diag_el) * sqrt(norm);
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
                double sum = 0;
                for (int k = it; k < N; ++k)
                {
                    sum += A[k*N + j] * v[k];
                }
                w[j] = beta * sum;
            }
            // Преобразование текущего блока матрицы A
            for (int j = it; j < N; ++j)
            {
                double* Aj = A + j * N;
                const double vj = v[j];
                for (int k = it; k < t; ++k)
                {
                    Aj[k] += vj * w[k];
                }
            }

            // * Заполнение поддиагнальной части it-го столбца матрицы A
            //   информативной частью вектора Хаусхолдера
            for (int j = it + 1; j < N; ++j)
            {
                A[j*N+it] = v[j];
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
                    Y[i*r] = v[i];
                    W[i*r] = beta * v[i];
                }
            }
            else
            {
                // Вычисление произведения (Y^t) * v
                for (int i = 0; i < it; ++i)
                {
                    double sum = 0;
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
                    double sum = 0;
                    for (int j = 0; j < it; ++j)
                    {
                        sum += Wi[j] * w[j];
                    }

                    buf = (i < shift) ? 0 : v[i];
                    W[i*r + it] = beta * (buf + sum);
                    Y[i*r + it] = buf;
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
            const double* Yi = Y + i*r;
            for (int j = lambda; j < N; ++j)
            {
                const double* Wj = W + j*r;
                double sum = 0;
                for (int k = 0; k < kbound; ++k)
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
    double* W = new double[N * b2];
    double* Y = new double[N * b2];
    double* WY = new double[N * N];
    double* Abuf = new double[N * N];

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
            double beta = A_diag_el + A_diag_el / abs(A_diag_el) * sqrt(norm);
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
                    w[j] += beta * sum;
                }
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
                        double sum = 0;
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
                    const int i_ub = min(mb_shift + b1, N) - mb_shift;
                    // высота блока ib
                    const int block_h = min(b1, N - mb_shift);
                    // указатели на начала блоков ib матриц Y и W
                    double* Yblock = Y + mb_shift*b2;
                    double* Wblock = W + mb_shift*b2;
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
            const int i_ub = min(i_mb_shift + b1, N) - i_mb_shift;
            const double* Yblock = Y + i_mb_shift*b2;

            const int ib_stripe_h = min(b1, N - i_mb_shift);
            double* WYstripe = WY + i_mb_shift*N + ib_stripe_h*d1_shift;

            for (int jb = 0; jb < col_b_count; ++jb)
            {
                const int j_mb_shift = d1_shift + jb*b1;
                const int j_lb = max(lambda, j_mb_shift) - j_mb_shift;
                const int j_ub = min(j_mb_shift + b1, N) - j_mb_shift;
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
                        double sum = 0;
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
            const int i_ub = min(i_mb_shift + b1, N) - i_mb_shift;
            const int ib_stripe_h = min(b1, N - i_mb_shift);

            // указатели сдвинуты до левой границы записываемой/читаемой области
            const double* Agen_ib_stripe_shifted = A + i_mb_shift*N + t*ib_stripe_h;
            double* Abuf_ib_stripe_shifted = Abuf + i_mb_shift*N + t*ib_stripe_h;
            const double* WY_ib_stripe_shifted = WY + i_mb_shift*N + d1_shift*ib_stripe_h;

            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_mb_shift = t + jb*b2;
                const int j_ub = min(j_mb_shift + b2, N) - j_mb_shift;
                const int jb_block_w = min(b2, N - j_mb_shift);

                // этот блок нужен для чтения на первой итерации по kb
                const double* Agen_init_ib_jb_block = Agen_ib_stripe_shifted + jb*ib_stripe_h*b2;
                // в этот блок будем записывать результат произведения
                double* Abuf_ib_jb_block = Abuf_ib_stripe_shifted + jb*ib_stripe_h*b2;

                // Для kb = 0 отдельно, чтобы не выполнять проверку
                // для присваивания в самом глубоком цикле (в А или в Abuf)
                const int k0_lb = max(lambda, d1_shift) - d1_shift;
                const int k0_ub = min(d1_shift + b1, N) - d1_shift;
                const int kb0_block_size = min(b1, N - d1_shift);

                const double* Aright_block_k0_jb = A + d1_shift*N + j_mb_shift*kb0_block_size;

                for (int i = i_lb; i < i_ub; ++i)
                {
                    double* Abuf_i = Abuf_ib_jb_block + i*jb_block_w;
                    const double* Agen_init_i = Agen_init_ib_jb_block + i*jb_block_w;
                    const double* WYi = WY_ib_stripe_shifted + i*kb0_block_size;

                    for (int j = 0; j < j_ub; ++j)
                    {
                        double sum = Agen_init_i[j];
                        for (int k = k0_lb; k < k0_ub; ++k)
                        {
                            sum += WYi[k] * Aright_block_k0_jb[k*jb_block_w + j];
                        }
                        Abuf_i[j] = sum;
                    }
                }

                for (int kb = 1; kb < col_b_count; ++kb)
                {
                    const int k_mb_shift = d1_shift + kb*b1;
                    const int kb_block_size = min(b1, N - k_mb_shift);
                    const int k_ub = min(k_mb_shift + b1, N) - k_mb_shift;
                    
                    const double* WYblock = WY_ib_stripe_shifted + kb*ib_stripe_h*b1;
                    const double* Ablock = A + k_mb_shift*N + j_mb_shift*kb_block_size;
                    for (int i = i_lb; i < i_ub; ++i)
                    {
                        double* Abuf_i = Abuf_ib_jb_block + i*jb_block_w;
                        const double* WYi = WYblock + i*kb_block_size;
                        for (int j = 0; j < j_ub; ++j)
                        {
                            double sum = 0;
                            for (int k = 0; k < k_ub; ++k)
                            {
                                sum += WYi[k] * Ablock[k*jb_block_w + j];
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
        const int ib_stripe_h = min(b1, N - d1_shift);

        const double* Abuf_stripe = Abuf + d1_shift*N;
        double* Astripe = A + d1_shift*N;
        for (int jb = 0; jb < row_b_count; ++jb)
        {
            const int j_mb_shift = t + jb*b2;
            const int jb_block_w = min(b2, N - j_mb_shift);

            double* A_block_shifted = Astripe + ib_stripe_h*j_mb_shift + i_lb*jb_block_w;
            const double* Abuf_block_shifted = Abuf_stripe + ib_stripe_h*j_mb_shift + i_lb*jb_block_w;
         
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
