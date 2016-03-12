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
                        double sum = 0;
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
