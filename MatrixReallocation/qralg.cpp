#include "qralg.h"
#include "service.h"

#include <algorithm>

using std::min;
using std::max;
using std::swap;

using std::cout;
using std::endl;////////////////////////////////////////////////////////////

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
    
    // lambda ��������� �� ������ �������� �����
    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // ��������� �� ������ ���������� ����� �� �������
        const int t = min(lambda + b2, N);
        const int alloc_size = (N - lambda)*sizeof(double);
        const int row_b_count = static_cast<int>(ceil(1.0 * (N - t) / b2));
        const int col_b_count = static_cast<int>(ceil(1.0 * (N - lambda) / b1));

        // ���������� �������������� ��� ������� ������
        for (int it = lambda; it < t; ++it)
        {
            // ������, ���� ��������� ��������� ������ ��������
            memset(w + lambda, 0, alloc_size);

            double norm = 0.0;
            double scalar = 1.0;
            double buf;
            // * ���������� ������� �����������
            for (int j = it; j < N; ++j)
            {
                buf = A[j * N + it];
                v[j] = buf;
                norm += buf * buf;
            }
            const double A_diag_el = A[it * N + it];
            double beta = A_diag_el + A_diag_el / abs(A_diag_el) * sqrt(norm);
            // ������������ ����������, �.�. v[it] = 1
            for (int j = it + 1; j < N; ++j)
            {
                buf = (v[j] /= beta);
                scalar += buf * buf;
            }
            v[it] = 1.0;

            // * ���������� ��������������
            beta = -2.0 / scalar;
            // ���������� ������� w = beta * (A(:,lambda:t-1)^t) * v .
            // �������� � it, �.�. it-�� ������� ����������� ��������� �� ������� A(it:,it:)
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
     
            // �������������� �������� ����� ������� A
            // ����� ��������� �� �����, ��� ���� �� ���������
            for (int j = it; j < N; ++j)
            {
                double* Aj = A + j * N;
                const double vj = v[j];
                for (int k = it; k < t; ++k)
                {
                    Aj[k] += vj * w[k];
                }
            }

            // * ���������� �������������� ����� it-�� ������� ������� A
            //   ������������� ������ ������� �����������
            for (int j = it + 1; j < N; ++j)
            {
                A[j*N + it] = v[j];
            }
        }// FOR

         // * ���������� WY-������������� ������������ ������ �����������
        const int d = t - lambda;
        for (int it = 0; it < d; ++it)
        {
            const int shift = lambda + it;
            double scalar = 1, buf;
            // * ���������� ����� � ��. ������������ ������� �����������
            for (int j = shift + 1; j < N; ++j)
            {
                buf = A[j * N + shift];
                scalar += buf * buf;
                v[j] = buf;
            }
            v[shift] = 1.0;
            double beta = -2.0 / scalar;

            // it - ����� ����������� �������� � �������� W � Y (= ������ ��������������� �������)
            if (it == 0)
            {
                // ��������� ���������� ������� �������
                // ���������: ������ lambda ����� ������ W � Y - �������,
                //            ������� ����� ����������� ���.
                for (int i = lambda; i < N; ++i)
                {
                    Y[i*b2] = v[i];
                    W[i*b2] = beta * v[i];
                }
            }
            else
            {
                memset(w, 0, d*sizeof(double));

                // ���������� ������������ (Y^t) * v
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

                // ���������� ������������ W * ((Y^t) * v)
                // ����� ���� ��������� �� ���������
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

        // * �������������� ��������� ����� ������� A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        // ��������� Y * W^t (������� WY ���������� ��������)
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

        // ��������� A(lamda:M-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_lb = lambda + ib*b1;
            const int i_ub = min(i_lb + b1, N);
            
            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_lb = t + jb*b2;
                const int j_ub = min(j_lb + b2, N);

                // ��� k = 0 ��������, ����� �� ��������� ��������
                // � ����� �������� �����
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
                
        // ����������� ��������������� ����� ������� A
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
    return NULL;
}

double* QR_WY_standard(double* A, const int N, const int r)
{
    double* v    = new double[N];
    double* w    = new double[N];
    double* W    = new double[N * r];
    double* Y    = new double[N * r];
    double* WY   = new double[N * N];
    double* Abuf = new double[N * N];

    // lambda ��������� �� ������ �������� �����
    for (int lambda = 0; lambda < N; lambda += r)
    {
        // ��������� �� ������ ���������� ����� �� �������
        const int t = min(lambda + r, N);

        // ���������� �������������� ��� ������� ������
        for (int it = lambda; it < t; ++it)
        {
            double norm = 0;
            double scalar = 1;
            double buf;
            // * ���������� ������� �����������
            for (int j = it; j < N; ++j)
            {
                buf = A[j * N + it];
                v[j] = buf;
                norm += buf * buf;
            }
            const double A_diag_el = A[it * N + it];
            double beta = A_diag_el + A_diag_el / abs(A_diag_el) * sqrt(norm);
            // ������������ ����������, �.�. v[it] = 1
            for (int j = it + 1; j < N; ++j)
            {
                buf = (v[j] /= beta);
                scalar += buf * buf;
            }
            v[it] = 1.0;

            // * ���������� ��������������
            beta = -2.0 / scalar;
            // ���������� ������� w = beta * (A(:,lambda:t-1)^t) * v .
            // �������� � it, �.�. it-�� ������� ����������� ��������� �� ������� A(it:,it:)
            for (int j = it; j < t; ++j)    
            {
                double sum = 0;
                for (int k = it; k < N; ++k)
                {
                    sum += A[k*N + j] * v[k];
                }
                w[j] = beta * sum;
            }
            // �������������� �������� ����� ������� A
            for (int j = it; j < N; ++j)
            {
                double* Aj = A + j * N;
                const double vj = v[j];
                for (int k = it; k < t; ++k)
                {
                    Aj[k] += vj * w[k];
                }
            }

            // * ���������� �������������� ����� it-�� ������� ������� A
            //   ������������� ������ ������� �����������
            for (int j = it + 1; j < N; ++j)
            {
                A[j*N+it] = v[j];
            }
        }// FOR
        
        // * ���������� WY-������������� ������������ ������ �����������
        const int d = t - lambda;
        for (int it = 0; it < d; ++it)
        {
            const int shift = lambda + it;
            double scalar = 1, buf;
            // * ���������� ����� � ��. ������������ ������� �����������
            for (int j = shift + 1; j < N; ++j)
            {
                buf = A[j * N + shift];
                scalar += buf * buf;
                v[j] = buf;
            }
            v[shift] = 1.0;
            double beta = -2.0 / scalar;

            // it - ����� ����������� �������� � �������� W � Y (= ������ ��������������� �������)
            if (it == 0)
            {
                // ��������� ���������� ������� �������
                // ���������: ������ lambda ����� ������ W � Y - �������,
                //            ������� ����� ����������� ���.
                for (int i = lambda; i < N; ++i)
                {
                    Y[i*r] = v[i];
                    W[i*r] = beta * v[i];
                }
            }
            else
            {
                // ���������� ������������ (Y^t) * v
                for (int i = 0; i < it; ++i)
                {
                    double sum = 0;
                    for (int j = shift; j < N; ++j)
                    {
                        sum += Y[j*r + i] * v[j];
                    }
                    w[i] = sum;
                }

                // ���������� ������������ W * ((Y^t) * v)
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

        // * �������������� ��������� ����� ������� A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)
        
        // ��������� Y * W^t (������� WY ���������� ��������)
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

        // ��������� A(lamda:M-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1)
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

        // ����������� ��������������� ����� ������� A
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
    // lambda ��������� �� ������ �������� �����
    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // ���������� �����, � ������� ��������� ������� (lambda,lambda).
        // ���������: i-���������� ����� ���������� ��-�� ������ �������� b1 � b2 ������. 
        curr_i_block_ind = lambda / b1;
        ++curr_j_block_ind;

        // ��������� �� ������ ���������� �� ����������� �����
        const int t = min(lambda + b2, N);
        // ���������� ����� ������ ������ �� �������� �� ������
        const int row_b_count = block_count_in_row - curr_j_block_ind - 1;
        // ���������� ����� ������ ������� ������� � ���� �� �������
        const int col_b_count = block_count_in_col - curr_i_block_ind;
        // ������ �� ��������� � ����������� ��� �������� �����
        const int d1_shift = curr_i_block_ind*b1;
        const int d2_shift = curr_j_block_ind*b2;

        // �������������� �������� ������/������ ������� ������� ������/�������
        const int d1 = min(N, (curr_i_block_ind + 1)*b1) - d1_shift;
        const int d2 = t - lambda;
        
        // * ���������� �������������� ��� ������� ������� ��������
        for (int it = 0; it < d2; ++it)
        {
            // ����� �� 'lambda+it'-�� ������� � �������
            const int full_shift = lambda + it;
            // �������� �������� � ��������� �� ��������� ����
            const int b_begin = (full_shift - d1_shift) / b1;

            memset(w, 0, d2*sizeof(double));
            memset(v + lambda, 0, it*sizeof(double));

            double norm = 0;
            double scalar = 1;
            double buf;

            // * ���������� ������� �����������
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                // ����� ����� ����� �� ����� �������
                const int mb_shift = jb*b1 + d1_shift;
                // ������� � ������ ������� ������ ����� � ����� jb (��������� ���������)
                const int j_lb = max(full_shift, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift+b1, N) - mb_shift;
                // ������ ����� jb
                const int block_h = min(b1, N - mb_shift);
                // ��������� �� ������ ����� jb
                const double* Ablock = A + mb_shift*N + block_h*d2_shift;
                double* v_jb_shift = v + mb_shift;
                // j - ����� ������ � ����� jb (��������� ���������)
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
            // ������������ ����������, �.�. v[it] = 1
            for (int j = full_shift + 1; j < N; ++j)
            {
                buf = (v[j] /= beta);
                scalar += buf * buf;
            }
            v[full_shift] = 1.0;
            beta = -2.0 / scalar;
        
            // * ���������� ��������������
            // ���������� ������� w = beta * (A(:,lambda:t-1)^t) * v .
            // �������� � it, �.�. it-�� ������� ����������� ��������� �� ������� A(it:,it:)
            for (int kb = b_begin; kb < col_b_count; ++kb)
            {
                // ����� ����� ����� �� ����� �������
                const int mb_shift = kb*b1 + d1_shift;
                // ������� � ������ ������� ������ ����� � ����� jb (��������� ���������)
                const int k_lb = max(full_shift, mb_shift) - mb_shift;
                const int k_ub = min(mb_shift + b1, N) - mb_shift;
                // ������ ����� jb
                const int block_h = min(b1, N - mb_shift);
                // ��������� �� ������ ����� jb
                const double* Ablock = A + mb_shift*N + block_h*d2_shift;
                const double* v_kb_shift = v + mb_shift;
                // j - ����� ������� � ����� kb (��������� ���������)
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

            // �������������� �������� ����� ������� A
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                // ����� ����� ����� �� ����� �������
                const int mb_shift = d1_shift + jb*b1;
                // ������� � ������ ������� ������ ����� � ����� jb (��������� ���������)
                const int j_lb = max(full_shift, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // ������ ����� jb
                const int j_block_h = min(b1, N - mb_shift);
                // ��������� �� ������ ����� jb
                double* Ablock = A + mb_shift*N + j_block_h*d2_shift;
                // j - ����� ������ � ����� jb (��������� ���������)
                for (int j = j_lb; j < j_ub; ++j)
                {
                    // ��������� �� ������ j-�� ������ ����� jb
                    double* Aj = Ablock + j*d2;
                    const double vj = v[mb_shift + j];
                    // k - ����� ������� � ����� jb (��������� ���������)
                    for (int k = it; k < d2; ++k)
                    {
                        Aj[k] += vj * w[k];
                    }
                }
            }
            
            // * ���������� �������������� ����� it-�� ������� ������� A
            //   ������������� ������ ������� �����������
            const int b_inc_begin = (full_shift+1 - curr_i_block_ind*b1) / b1;
            for (int jb = b_inc_begin; jb < col_b_count; ++jb)
            {
                // ����� ����� ����� �� ����� �������
                const int mb_shift = d1_shift + jb*b1;
                // ������� � ������ ������� ������ ����� � ����� jb (��������� ���������)
                const int j_lb = max(full_shift+1, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // ������ ����� jb
                const int j_block_h = min(b1, N - mb_shift);
                // ��������� �� ������ ����� jb
                double* Ablock = A + mb_shift*N + j_block_h*d2_shift;
                const double* v_shift = v + mb_shift;
                // j - ����� ������ � ����� jb (��������� ���������)
                for (int j = j_lb; j < j_ub; ++j)
                {
                    Ablock[j*d2 + it] = v_shift[j];
                }
            }
        }// FOR


        memset(W, 0, N*b2*sizeof(double));
        memset(Y, 0, N*b2*sizeof(double));

        // * ���������� WY-������������� ������������ ������ �����������
        for (int it = 0; it < d2; ++it)
        {
            // ����� �� 'lambda+it'-�� ������� � �������
            const int full_shift = lambda + it;
            // �������� �������� � ��������� �� ��������� ����
            const int b_begin = (full_shift - d1_shift) / b1;

            memset(w, 0, d2*sizeof(double));
            memset(v + lambda, 0, it*sizeof(double));

            double scalar = 1, buf;
            // * ���������� ����� � ��. ������������ ������� �����������
            for (int jb = (full_shift + 1 - d1_shift) / b1; jb < col_b_count; ++jb)
            {
                // ����� ����� ����� �� ����� �������
                const int mb_shift = jb*b1 + d1_shift;
                // ������� � ������ ������� ������ ����� � ����� jb (��������� ���������)
                const int j_lb = max(full_shift+1, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // ������ ����� jb
                const int block_h = min(b1, N - mb_shift);
                // ��������� �� ������ ����� jb
                const double* Ablock = A + mb_shift*N + block_h*d2_shift;
                double* v_jb_shift = v + mb_shift;
                // j - ����� ������ � ����� jb (��������� ���������)
                for (int j = j_lb; j < j_ub; ++j)
                {
                    buf = Ablock[j*d2 + it];
                    scalar += buf * buf;
                    v_jb_shift[j] = buf;
                }
            }
            v[full_shift] = 1.0;
            double beta = -2.0 / scalar;
            
            // it - ����� ����������� �������� � �������� W � Y (= ������ ��������������� �������)
            if (it == 0)
            {
                // ��������� ���������� ������� �������
                // ���������: ������ lambda ����� ������ W � Y - �������,
                //            ������� ����� ����������� ���.
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    // ����� ����� ����� �� ����� �������
                    const int mb_shift = ib*b1 + d1_shift;
                    // ������� � ������ ������� ������ ����� � ����� ib (��������� ���������)
                    const int i_lb = max(full_shift, mb_shift) - mb_shift;
                    const int i_ub = min(mb_shift + b1, N) - mb_shift;
                    // ������ ����� ib
                    const int block_h = min(b1, N - mb_shift);
                    // ��������� �� ������ ������ ib ������ Y � W
                    double* Yblock = Y + mb_shift*b2;
                    double* Wblock = W + mb_shift*b2;
                    const double* v_shift = v + mb_shift;
                    for (int i = i_lb; i < i_ub; ++i)
                    {
                        Yblock[i*b2] = v_shift[i];
                        Wblock[i*b2] = beta * v_shift[i];
                    }
                }
//                print_to(cout, Y, N, b2);
//                print_to(cout, W, N, b2);
            }
            else
            {
                // ���������� ������������ (Y^t) * v
                for (int jb = b_begin; jb < col_b_count; ++jb)
                {
                    // ����� ����� ����� �� ����� �������
                    const int mb_shift = jb*b1 + d1_shift;
                    // ������� � ������ ������� ������ ����� � ����� jb (��������� ���������)
                    const int j_lb = max(full_shift, mb_shift) - mb_shift;
                    const int j_ub = min(mb_shift + b1, N) - mb_shift;
                    // ������ ����� jb
                    const int block_h = min(b1, N - mb_shift);
                    // ��������� �� ������ ����� jb
                    const double* Yblock = Y + mb_shift*b2;
                    double* v_shift = v + mb_shift;
                    // j - ����� ������ � ����� jb (��������� ���������)
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
                
                // ���������� ������������ W * ((Y^t) * v)
                // ����� ���� ��������� �� ���������
                for (int ib = 0; ib < col_b_count; ++ib)
                {
                    // ����� ����� ����� �� ����� �������
                    const int mb_shift = ib*b1 + d1_shift;
                    // ������� � ������ ������� ������ ����� � ����� ib (��������� ���������)
                    const int i_lb = max(lambda, mb_shift) - mb_shift;
                    const int i_ub = min(mb_shift + b1, N) - mb_shift;
                    // ������ ����� ib
                    const int block_h = min(b1, N - mb_shift);
                    // ��������� �� ������ ������ ib ������ Y � W
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

  //      memset(WY, 0, N*N*sizeof(double));
        
        // * �������������� ��������� ����� ������� A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)
        
        // ��������� Y * W^t (������� WY ���������� ��������).
        // ������� WY ����� ����� b1 x b1 ����������
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
                    // �����������, ��������� �� ����������� ����� ������� Y
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

        // ��������� A(lamda:M-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1)
        // ����������: (b1 x b2) = (b1 x b1) * (b1 x b2)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int i_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int i_ub = min(i_mb_shift + b1, N) - i_mb_shift;
            const int ib_stripe_h = min(b1, N - i_mb_shift);

            // ��������� �������� �� ����� ������� ������������/�������� �������
            const double* Agen_ib_stripe_shifted = A + i_mb_shift*N + t*ib_stripe_h;
            double* Abuf_ib_stripe_shifted = Abuf + i_mb_shift*N + t*ib_stripe_h;
            const double* WY_ib_stripe_shifted = WY + i_mb_shift*N + d1_shift*ib_stripe_h;

            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_mb_shift = t + jb*b2;
                const int j_ub = min(j_mb_shift + b2, N) - j_mb_shift;
                const int jb_block_w = min(b2, N - j_mb_shift);

                // ���� ���� ����� ��� ������ �� ������ �������� �� kb
                const double* Agen_init_ib_jb_block = Agen_ib_stripe_shifted + jb*ib_stripe_h*b2;
                // � ���� ���� ����� ���������� ��������� ������������
                double* Abuf_ib_jb_block = Abuf_ib_stripe_shifted + jb*ib_stripe_h*b2;

                // ��� kb = 0 ��������, ����� �� ��������� ��������
                // ��� ������������ � ����� �������� ����� (� � ��� � Abuf)
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
        
        // * ����������� ��������������� ����� ������� A
        // ����������� ������ ������� ������ (��������, �.�. ��� ����� ���� ������� �� ������)
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
    double* sum_vec = new double[b1];

    int curr_i_block_ind = -1;
    int curr_j_block_ind = -1;

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
            // ������ �������� �����, ����������� 'lambda+it'-�� ������������ �������,
            // ������ �� �����, ����������� 'lambda'-�� �������
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
                    // Shift of current small-block-row from matrix first row
                    const int bd_shift = id*db1 + mb_shift;
                    // Bounds for current small block reading.
                    // Reasons of using are the same as 'i_lb' and 'i_ub'
                    const int lb = max(full_shift, bd_shift) - bd_shift;
                    const int ub = min(bd_shift + db1, i_ub) - bd_shift;
                    // 'id'-th small block height (number of rows)
                    const int id_block_h = min(db1, ib_block_h - id*db1);
                    // Pointer to beginning of required small block in 'id'-th stripe
                    const double* Adblock = Ablock + id*db1*d2 +
                        dblock_col_index*id_block_h*db2;
                    double* v_shifted = v + bd_shift;
                    
                    for (int i = lb; i < ub; ++i)
                    {
                        buf = Adblock[i*dblock_width + loc_it];
                        v_shifted[i] = buf;
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

            double beta = A_diag_el + A_diag_el / abs(A_diag_el) * sqrt(norm);
            // Use normalization so, that v[it] = 1
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
                
                const int kd_lb = k_lb / db1;
                const int kd_ub = (kb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int kd = kd_lb; kd < kd_ub; ++kd)
                {
                    const int bd_shift = kd*db1 + mb_shift;
                    const int lb = max(full_shift, bd_shift) - bd_shift;
                    const int ub = min(bd_shift + db1, k_ub) - bd_shift;
                    const int dblock_h = min(db1, block_h - kd*db1);
                    // Pointer to beginning of 'kb'-th small block stripe
                    const double* Adstripe = Ablock + kd*db1*d2;
                    const double* v_kd_shift = v_kb_shift + kd*db1;
                    
                    // const int& jd_lb = dblock_col_index;
                    // const int& jd_ub = dblock_count_in_curr_block_row;
                    for (int jd = dblock_col_index; jd < dblock_count_in_curr_block_row; ++jd)
                    {
                        // Width of this small block
                        const int dblock_w = min(db2, d2 - jd*db2);
                        const int loc_drow_shift = jd * db2;
                        const int j_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int j_ub = min(d2, loc_drow_shift + db2) - loc_drow_shift;
                        // Pointer to beginning of 'jb'-th small block
                        const double* Adblock = Adstripe + jd*dblock_h*db2;
                        double* w_shifted = w + loc_drow_shift;

                        for (int j = j_lb; j < j_ub; ++j)
                        {
                            double sum = 0.0;
                            for (int k = lb; k < ub; ++k)
                            {
                                sum += Adblock[k*dblock_w + j] * v_kd_shift[k];
                            }
                            w_shifted[j] += beta * sum;
                        }
                    }
                }
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

                const int jd_lb = jb_lb / db1;
                const int jd_ub = (jb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int jd = jd_lb; jd < jd_ub; ++jd)
                {
                    const int bd_shift = jd*db1 + mb_shift;
                    const int j_lb = max(full_shift, bd_shift) - bd_shift;
                    const int j_ub = min(bd_shift + db1, jb_ub) - bd_shift;
                    const int jd_block_h = min(db1, jb_block_h - jd*db1);
                    double* Adstripe = Ablock + jd*db1*d2;
                    const double* v_shifted = v + bd_shift;

                    for (int kd = dblock_col_index; kd < dblock_count_in_curr_block_row; ++kd)
                    {
                        const int loc_drow_shift = kd * db2;
                        const int k_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int k_ub = min(d2, loc_drow_shift + db2) - loc_drow_shift;
                        const int dblock_w = min(db2, d2 - kd*db2);
                        double* Adblock = Adstripe + kd*jd_block_h*db2;
                        const double* w_shifted = w + loc_drow_shift;

                        for (int j = j_lb; j < j_ub; ++j)
                        {
                            double* Aj = Adblock + j*dblock_w;
                            const double vj = v_shifted[j];
                            for (int k = k_lb; k < k_ub; ++k)
                            {
                                Aj[k] += vj * w_shifted[k];
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
                    const int bd_shift = id*db1 + mb_shift;
                    const int i_lb = max(full_shift + 1, bd_shift) - bd_shift;
                    const int i_ub = min(bd_shift + db1, ib_ub) - bd_shift;
                    const int id_block_h = min(db1, ib_block_h - id*db1);
                    double* Adblock = Ablock + id*db1*d2 +
                        dblock_col_index * id_block_h * db2;
                    const double* v_shifted = v + bd_shift;
                    
                    for (int i = i_lb; i < i_ub; ++i)
                    {
                        Adblock[i*dblock_width + loc_it] = v_shifted[i];
                    }
                }
            }
        }// FOR
        
        memset(W, 0, N*b2*sizeof(double));
        memset(Y, 0, N*b2*sizeof(double));

        // * Computation of WY-representation of Householder matrices production...
//        memset(W + d1_shift*b2, 0, min(bsizes_ratio*b1, N - d1_shift)*b2*sizeof(double));
//        memset(Y + d1_shift*b2, 0, min(bsizes_ratio*b1, N - d1_shift)*b2*sizeof(double));
        for (int it = 0; it < d2; ++it)
        {
            memset(w, 0, d2*sizeof(double));
            memset(v + lambda, 0, it*sizeof(double));

            const int full_shift = lambda + it;
            const int b_begin = (full_shift - d1_shift) / b1;
            const int dblock_col_index = it / db2;
            const int dblock_row_index = (full_shift - (d1_shift + b_begin*b1)) / db1;
            const int dblock_width = min(db2, d2 - dblock_col_index*db2);
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
                    const int bd_shift = mb_shift + id*db1;
                    const int lb = max(full_shift + 1, bd_shift) - bd_shift;
                    const int ub = min(bd_shift + db1, i_ub) - bd_shift;
                    const int id_block_h = min(db1, ib_block_h - id*db1);
                    const double* Adblock = Ablock + id*db1*d2 +
                        dblock_col_index*id_block_h*db2;
                    double* v_shifted = v + bd_shift;

                    for (int i = lb; i < ub; ++i)
                    {
                        buf = Adblock[i*dblock_width + loc_it];
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
                    const int ib_block_h = min(b1, N - mb_shift);

                    const int id_lb = (i_lb == 0) ? 0 : (i_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;
                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int bd_shift = id*db1 + mb_shift;
                        const int lb = max(lambda, bd_shift) - bd_shift;
                        const int ub = min(bd_shift + db1, i_ub) - bd_shift;
                        const int id_block_h = min(db1, ib_block_h - id*db1);
                        // Width of block is 'b2',
                        // because memory was allocated for complete matrices
                        double* Ydblock = Y + bd_shift*b2;
                        double* Wdblock = W + (Ydblock - Y);
                        double* v_shifted = v + bd_shift;

                        for (int i = lb; i < ub; ++i)
                        {
                            Ydblock[i*dblock_width_WY] = v_shifted[i];
                            Wdblock[i*dblock_width_WY] = beta * v_shifted[i];
                        }
                    }
                }
 //               print_to(cout, Y, N, b2);
 //               print_to(cout, W, N, b2);
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
                        const int bd_shift = id*db1 + mb_shift;
                        const int i_lb = max(full_shift, bd_shift) - bd_shift;
                        const int i_ub = min(bd_shift + db1, ib_ub) - bd_shift;
                        const int id_block_h = min(db1, ib_block_h - id*db1);
                        const double* Ydstripe = Y + bd_shift*b2;
                        const double* v_shifted = v + bd_shift;

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it, loc_drow_shift + db2) - loc_drow_shift;
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

                memset(WY, 0, N*N*sizeof(double));

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
                        const int bd_shift = id*db1 + mb_shift;
                        const int i_lb = max(lambda, bd_shift) - bd_shift;
                        const int i_ub = min(bd_shift + db1, ib_ub) - bd_shift;
                        const int id_block_h = min(db1, ib_block_h - id*db1);
                        double* Wdstripe = W + bd_shift*b2;
                        double* Ydstripe = Y + (Wdstripe - W);
                        // Vector to contain results of production of W's row and vector w
                        memset(sum_vec, 0, b1*sizeof(double));

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
                        const double* v_shifted = v + bd_shift;
                        // Pointers to small blocks containing 'it'-th column of Y and W
                        double* Wdblock_add = Wdstripe + dblock_col_index*id_block_h*db2;
                        double* Ydblock_add = Ydstripe + (Wdblock_add - Wdstripe);
                        for (int i = i_lb; i < i_ub; ++i)
                        {
                            buf = v_shifted[i];
                            Wdblock_add[i*dblock_width_WY + loc_it] = beta * (buf + sum_vec[i]);
                            Ydblock_add[i*dblock_width_WY + loc_it] = buf;
                        }
                    }
                }
            }
        }

        // * Transformation of remaining block-columns of A...
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        const int up_stripe_h = min(b1, N - d1_shift);
        memset( WY + d1_shift*N + d1_shift*up_stripe_h,
                0,
                ((N-d1_shift)*N - d1_shift*up_stripe_h)*sizeof(double));

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
                    const int i_bd_shift = i_mb_shift + id*db1;
                    const int i_lb = max(lambda, i_bd_shift) - i_bd_shift;
                    const int i_ub = min(i_bd_shift + db1, ib_ub) - i_bd_shift;
                    const int id_block_h = min(db1, ib_stripe_h - id*db1);
                    const double* Ydstripe = Y + i_bd_shift*b2;
                    double* WYdstripe = WYblock + id*db1*jb_stripe_h;
                    
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
                        const int j_bd_shift = j_mb_shift + jd*db1;
                        const int j_lb = max(lambda, j_bd_shift) - j_bd_shift;
                        const int j_ub = min(j_bd_shift + db1, jb_ub) - j_bd_shift;
                        const int jd_block_h = min(db1, jb_stripe_h - jd*db1);
                        const double* Wdstripe = W + j_bd_shift*b2;
                        double* WYdblock = WYdstripe + jd*id_block_h*db1;
                        
                        for (int kd = 0; kd < kd_ud; ++kd)
                        {
                            const int loc_drow_shift = kd * db2;
                            const int k_ub = min(d2, loc_drow_shift + db2) - loc_drow_shift;
                            const int kd_dblock_w = min(db2, b2 - loc_drow_shift);
                            const double* Ydblock = Ydstripe + id_block_h*loc_drow_shift;
                            const double* Wdblock = Wdstripe + jd_block_h*loc_drow_shift;

                            for (int i = i_lb; i < i_ub; ++i)
                            {
                                double* WYrow = WYdblock + i*jd_block_h;
                                const double* Yrow = Ydblock + i*kd_dblock_w;
                                
                                for (int j = j_lb; j < j_ub; ++j)
                                {
                                    const double* Wrow = Wdblock + j*kd_dblock_w;
                                    double sum = 0.0;
                                    for (int k = 0; k < k_ub; ++k)
                                    {
                                        sum += Yrow[k] * Wrow[k];
                                    }
                                    WYrow[j] += sum;
                                }
                            }
                        }
                    }
                }
            }
        }

        memcpy( Abuf + d1_shift*N + d2_shift*up_stripe_h,
                A + d1_shift*N + d2_shift*up_stripe_h,
                ((N - d1_shift)*N - d2_shift*up_stripe_h)*sizeof(double));

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
            double* Abuf_stripe = Abuf + i_mb_shift*N;

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
                        const int i_bd_shift = i_mb_shift + id*db1;
                        const int i_lb = max(lambda, i_bd_shift) - i_bd_shift;
                        const int i_ub = min(i_bd_shift + db1, ib_ub) - i_bd_shift;
                        const int id_block_h = min(db1, ib_stripe_h - id*db1);
                        const double* WY_dstripe = WY_block + id*db1*kb_stripe_h;
                        double* Abuf_dstripe = Abuf_block + id*db1*jb_column_w;

                        for (int jd = 0; jd < jd_ub; ++jd)
                        {
                            const int j_bd_shift = j_mb_shift + jd*db2;
                            const int j_ub = min(j_bd_shift + db2, jb_ub) - j_bd_shift;
                            const int jd_column_w = min(db2, jb_column_w - jd*db2);
                            double* Abuf_dblock = Abuf_dstripe + jd*id_block_h*db2;

                            // Remaining 'kd'-loop iterations
                            for (int kd = kd_lb; kd < kd_ub; ++kd)
                            {
                                const int k_bd_shift = k_mb_shift + kd*db1;
                                const int k_ub = min(k_bd_shift + db1, kb_ub) - k_bd_shift;
                                const int kd_block_h = min(db1, kb_stripe_h - kd*db1);
                                const double* WY_dblock = WY_dstripe + kd*id_block_h*db1;
                                const double* A_dblock = A_block + kd*db1*jb_column_w + jd*kd_block_h*db2;

                                for (int i = i_lb; i < i_ub; ++i)
                                {
                                    double* Abuf_row = Abuf_dblock + i*jd_column_w;
                                    const double* WY_row = WY_dblock + i*kd_block_h;
                                    for (int j = 0; j < j_ub; ++j)
                                    {
                                        double sum = 0.0;
                                        for (int k = 0; k < k_ub; ++k)
                                        {
                                            sum += WY_row[k] * A_dblock[k*jd_column_w + j];
                                        }
                                        Abuf_row[j] += sum;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Copying of transformed part of A matrix, that is stored in Abuf matrix...
        // Firstly, copy of upper big-block-stripe, because it can be truncated from above
        const int ib_lb = lambda - d1_shift;                // local
        const int ib_ub = min(d1_shift + b1, N) - d1_shift; // local
        const int id_border_index = ib_lb / db1;
        const int i_bd_shift = d1_shift + id_border_index * db1;
        const int i_lb = lambda - i_bd_shift;
        const int i_ub = min(i_bd_shift + db1, N) - i_bd_shift;
        const int stripe_h = min(b1, N - d1_shift);
        const int dstripe_h = min(db1, stripe_h - id_border_index*db1);

        double* A_stripe = A + d1_shift*N;
        const double* Abuf_stripe = Abuf + d1_shift*N;
        for (int jb = 0; jb < row_b_count; ++jb)
        {
            const int j_mb_shift = t + jb*b2;
            const int jb_block_w = min(b2, N - j_mb_shift);

            const double* Abuf_block = Abuf_stripe + j_mb_shift*stripe_h;
            const double* Abuf_dstripe = Abuf_block + id_border_index*db1*jb_block_w;
            double* A_block = A_stripe + (Abuf_block - Abuf_stripe);
            double* A_dstripe = A_block + (Abuf_dstripe - Abuf_block);
            
            const int jd_ub = (jb == row_b_count - 1) ?
                dblock_count_in_diff_brow :
                dblock_count_in_brow;
            for (int jd = 0; jd < jd_ub; ++jd)
            {
                const int j_bd_shift = j_mb_shift + jd*db2;
                const int jd_block_w = min(db2, jb_block_w - jd*db2);
                const double* Abuf_dblock_shifted = Abuf_dstripe +
                    jd*dstripe_h*db2 +      // shift to small block
                    i_lb*jd_block_w;    // shift to required row
                double* A_dblock_shifted = A_dstripe +
                    (Abuf_dblock_shifted - Abuf_dstripe);

                memcpy(A_dblock_shifted, Abuf_dblock_shifted, (i_ub - i_lb)*jd_block_w*sizeof(double));
            }

            // Copying of 'jb'-th big block remaining elements
            const int shift = (id_border_index + 1) * db1;
            if (shift < stripe_h)
            {
                memcpy(A_dstripe + shift * jb_block_w,
                       Abuf_dstripe + shift * jb_block_w,
                      (stripe_h - shift) * jb_block_w * sizeof(double));
            }
        }

        for (int ib = 1; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int ib_stripe_h = min(b1, N - i_mb_shift);
            // Pointers to first block of big-block-stripe, which must be copied
            const double* Abuf_stripe_shifted = Abuf + i_mb_shift*N + t*ib_stripe_h;
            double* A_stripe_shifted = A + (Abuf_stripe_shifted - Abuf);

            memcpy(A_stripe_shifted, Abuf_stripe_shifted, ib_stripe_h*(N - t)*sizeof(double));
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
