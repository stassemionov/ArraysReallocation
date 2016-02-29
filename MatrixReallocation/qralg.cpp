#include "qralg.h"

#include <algorithm>

using std::min;
using std::swap;

double* HouseholderWY(double* A, const int M, const int N, const int r)
{
    double* v = new double[M];
    double* w = new double[N];
    double* W = new double[M * r];
    double* Y = new double[M * r];
    double* BUF = new double[(N - r) * r];
    double* Yv = new double[r]; // Yv - это произведение матрицы transp(Y) и вектора w

    int lambda = 0;
    while (lambda < N - 1)
    {
        // lambda указывает на начало текущего блока, t - на начало следующего за ним
        int t = min(lambda + r, N - 1);

        // ** Выполнение преобразований над текущим блоком
        for (int it = lambda; it < t; ++it)
        {
            double norm = 0;
            double scalar = 1;
            double buf = 0;
            // * Вычисление вектора Хаусхолдера
            for (int j = it; j < M; ++j)
            {
                buf = A[j * N + it];
                v[j] = buf;
                norm += buf * buf;
            }
            buf = A[it * N + it];
            double beta = buf + buf / abs(buf) * sqrt(norm);
            // используется нормировка т.ч. v[it] = 1
            for (int j = it + 1; j < M; ++j)
            {
                buf = v[j] / beta;
                scalar += buf * buf;
                v[j] = buf;
            }
            v[it] = 1;

            // * Применение преобразования
            beta = -2 / scalar;
            // вычисление вектора w 
            for (int j = it; j < t; ++j)
            {
                double sum = 0;
                for (int k = it; k < M; ++k)
                {
                    sum += A[k*N + j] * v[k];
                }
                w[j] = beta * sum;
            }
            // преобразование текущего блока матрицы A
            for (int j = it; j < M; ++j)
            {
                double* AjN = A + j * N;
                for (int k = it; k < t; ++k)
                {
                    AjN[k] += v[j] * w[k];
                }
            }

            // * Заполнение поддиагнальной части столбца матрицы A информативной частью вектора Хаусхолдера
            for (int j = it + 1; j < M; ++j)
            {
                A[j * N + it] = v[j];
            }

            // * Вычисление WY-представления произведения матриц Хаусхолдера
            // l - число уже накопленных столбцов в матрицах W и Y (~ номеру дозаписываемого столбца)
            const int l = it - lambda;
            if (l == 0)
            {
                // начальное заполнение первого столбца
                for (int i = it; i < M; ++i)
                {
                    Y[i * r] = v[i];
                    W[i * r] = beta * v[i];
                }
            }
            else	// Замечание: первые lambda строк матриц W и Y - нулевые
            {
                // матрицы W и Y уже содержат l столбцов.
                // Умножение transp(Y) на вектор v
                for (int i = 0; i < l; ++i)
                {
                    double sum = 0;
                    for (int j = it; j < M; ++j)
                    {
                        sum += Y[j * r + i] * v[j]; /////////////////////////////////////
                    }
                    Yv[i] = sum;
                }

                //
                for (int i = lambda; i < M; ++i)
                {
                    double sum = 0;
                    double* Wir = W + i * r;
                    for (int j = 0; j < l; ++j)
                    {
                        sum += Wir[j] * Yv[j];
                    }

                    // Проверка нужна, чтобы не использовать неактуальные координаты вектора v
                    if (i >= it)
                    {
                        // дозапись столбцов в матрицы W и Y
                        Wir[l] = beta * (v[i] + sum);
                        Y[i*r + l] = v[i];
                    }
                    else
                    {
                        // Здесь v[i] равно 0, поэтому отсутствует
                        Wir[l] = beta * sum;
                    }
                }
            }
        }// FOR

        const int d = t - lambda;
        // ** произведение выполняется транспонированно ( форма: (N-r) x r )
        for (int i = t; i < N; ++i)	// строка transp(A) ~ столбцу A
        {
            double* bsh = BUF + (i - r) * r;
            for (int j = 0; j < d; ++j)	// столбец W
            {
                double sum = 0;
                for (int k = lambda; k < M; ++k)	// столбец transp(A) ~ строке A
                {
                    sum += A[k * N + i] * W[k * r + j];
                }
                bsh[j] = sum;
            }
        }

        for (int i = lambda; i < M; ++i)
        {
            double* Yir = Y + i*r;
            for (int j = t; j < N; ++j)
            {
                double* BUFjr = BUF + (j - r)*r;
                const int ub = min(i - lambda + 1, d);
                double sum = 0;
                for (int k = 0; k < ub; ++k)
                {
                    sum += Yir[k] * BUFjr[k];
                }
                A[i*N + j] += sum;
            }
        }

        lambda = t;
    }// WHILE

    return A;
}

double* get_inverse_matrix(double* A, int N)
{
    double* B = new double[N*N], *B_sh = 0;
    for (int i = 0; i < N; i++)
    {
        B_sh = B + i*N;
        for (int j = 0; j < N; j++)
        {
            B_sh[j] = (i == j);
        }
    }

    int* used = new int[N];	// used: #итерации -> #строки
    for (int i = 0; i < N; i++)
        used[i] = -1;

    bool was_used = false;
    int index = -1;
    double main_el, buf;
    double* A_i_sh, *A_k_sh, *B_i_sh, *B_k_sh;
    for (int it = 0; it < N; it++)	// номер итерации (номер столбца)
    {
        for (int i = 0; i < N; i++)	// выбор ведущей строки
        {
            was_used = false;
            for (int k = 0; k < it; k++)	// в любом случае, нужно проверить, выбиралась ли ранее эта строка ведущей 
            {
                if (used[k] == i)
                {
                    was_used = true;
                    break;
                }
            }
            if ((A[i*N + it] == 0) || was_used)	// если делящий элемент = 0, либо эта строка уже была ведущей
            {
                index = -1;
                for (int k = 0; k < N; k++)	// поиск ненулевого элемента по столбцу в неиспользованных строках
                {
                    if ((A[k*N + it] != 0) && (used[k] == -1))
                    {
                        index = k;
                        break;
                    }
                }
                if (index == -1) // весь столбец нулевой
                    continue;	// переход к следующему столбцу
            }
            else
            {
                index = i;
            }

            used[it] = index;
            main_el = A[index*N + it];

            B_i_sh = B + index*N;
            A_i_sh = A + index*N;
            for (int j = 0; j < it; j++)
            {
                B_i_sh[j] /= main_el;
            }
            for (int j = it; j < N; j++)
            {
                A_i_sh[j] /= main_el;
                B_i_sh[j] /= main_el;
            }
            for (int k = 0; k < N; k++)	// все строки
            {
                if (k != index)	// кроме ведущей
                {
                    B_k_sh = B + k*N;
                    A_k_sh = A + k*N;
                    buf = A[k*N + it];
                    for (int j = 0; j < it; j++)
                    {
                        B_k_sh[j] -= B_i_sh[j] * buf;
                    }
                    for (int j = it; j < N; j++)
                    {
                        A_k_sh[j] -= A_i_sh[j] * buf;
                        B_k_sh[j] -= B_i_sh[j] * buf;
                    }
                }
            }
        }
    }
    // восстановление порядка строк
    for (int i = 0; i < N; i++)
    {
        B_k_sh = B + used[i] * N;
        B_i_sh = B + i*N;
        for (int j = 0; j < N; j++)
        {
            swap(B_k_sh[j], B_i_sh[j]);
        }
    }
    return B;
}
