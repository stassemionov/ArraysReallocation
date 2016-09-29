#include "floydalg.h"

#include <algorithm>

using std::min;

double* floyd_standard(double* data, const int N)
{
    for (int k = 0; k < N; ++k)
    {
        const double* Ak = data + k * N;
        for (int i = 0; i < N; ++i)
        {
            double* Ai = data + i * N;
            const double Aik = Ai[k];
            for (int j = 0; j < N; ++j)
            {
                Ai[j] = min(Ai[j], Aik + Ak[j]);
            }
        }
    }
    return data;
}

double* floyd_tiled(double* data, const int N, const int B)
{
    const int factor = static_cast<int>(ceil(1.0 * N / B));

    // Шагаем по блокам (factor - их количество)
    for (int t = 0; t < factor; ++t)
    {
        const int tmB = t * B;
        const int tB = min(tmB + B, N);
        //* вычисления фазы 1 *//
        //* блок на перекрестье строки и столбца *//
        for (int k = tmB; k < tB; ++k)
        {
            const double* Ak = data + k * N;
            for (int i = tmB; i < tB; ++i)
            {
                double* Ai = data + i * N;
                const double Aik = Ai[k];
                for (int j = tmB; j < tB; ++j)
                {
                    Ai[j] = min(Ai[j], Aik + Ak[j]);
                }
            }
        }

        //* вычисления фазы 2 *//
        //* вычисления блоков по строке *//
        
        for (int b = 0; b < factor; ++b)
        {
            if (b == t)
                continue;

            const int lb = b * B;
            const int ub = min(lb + B, N);

            for (int k = tmB; k < tB; ++k)
            {
                const double* Ak = data + k * N;
                for (int i = tmB; i < tB; ++i)
                {
                    double* Ai = data + i * N;
                    double Aik = Ai[k];
                    for (int j = lb; j < ub; ++j)
                    {
                        Ai[j] = min(Ai[j], Aik + Ak[j]);
                    }
                }
            }
        }
        //* вычисления блоков по столбцу *//
        for (int b = 0; b < factor; ++b)
        {
            if (b == t)
                continue;

            const int lb = b * B;
            const int ub = min(lb + B, N);
            for (int k = tmB; k < tB; ++k)
            {
                const double* Ak = data + k * N;
                for (int i = lb; i < ub; ++i)
                {
                    double* Ai = data + i * N;
                    const double Aik = Ai[k];
                    for (int j = tmB; j < tB; ++j)
                    {
                        Ai[j] = min(Ai[j], Aik + Ak[j]);
                    }
                }
            }
        }

        //* вычисления фазы 3 *//
        for (int b1 = 0; b1 < factor; ++b1)
        {
            if (b1 == t)
                continue;

            const int lb1 = b1 * B;
            const int ub1 = min(lb1 + B, N);
            for (int b2 = 0; b2 < factor; ++b2)
            {
                if (b2 == t)
                    continue;

                const int lb2 = b2 * B;
                const int ub2 = min(lb2 + B, N);
                for (int k = tmB; k < tB; ++k)
                {
                    const double* Ak = data + k * N;
                    for (int i = lb1; i < ub1; ++i)
                    {
                        double* Ai = data + i * N;
                        const double Aik = Ai[k];
                        for (int j = lb2; j < ub2; ++j)
                        {
                            Ai[j] = min(Ai[j], Aik + Ak[j]);
                        }
                    }
                }
            }
        }
    }
    return data;
}

double* floyd_block(double* data, const int N, const int B)
{
    const int factor = static_cast<int>(ceil(1.0 * N / B));

    for (int t = 0; t < factor; ++t)
    {
        // Размер центрального блока (на пересечений t-го столбца и t-ой строки).
        // центральный блок всегда диагональный
        const int CBsize = min(B, N - B*t);

        // указатель на начало центрально блока
        double* CB = data + t*B*N + t*CBsize*B;

        for (int k = 0; k < CBsize; ++k)
        {
            const double* CBk = CB + k*CBsize;
            for (int i = 0; i < CBsize; ++i)
            {
                double* Ai = CB + i*CBsize;
                const double Aik = Ai[k];
                for (int j = 0; j < CBsize; ++j)
                {
                    Ai[j] = min(Ai[j], Aik + CBk[j]);
                }
            }
        }

        //* вычисления блоков по строке *//
        for (int b = 0; b < factor; ++b)
        {
            if (b == t)
                continue;

            // ширина текущего (b-го по t-ой строке) блока
            const int ABwidth = min(B, N - b*B);
            // указатель на начало текущего блока
            double* ABi = data + t*B*N + b*CBsize*B;
            for (int k = 0; k < CBsize; ++k)
            {
                // указатель на начало к-ой строки текущего блока
                const double* Ak = ABi + k*ABwidth;
                for (int i = 0; i < CBsize; ++i)
                {
                    // указатель на i-ую строку текущего блока
                    double* Ai = ABi + i*ABwidth;
                    // к-ый элемент i-ой строки центрального блока
                    const double Aik = CB[i*CBsize + k];
                    for (int j = 0; j < ABwidth; ++j)
                    {
                        Ai[j] = min(Ai[j], Aik + Ak[j]);
                    }
                }
            }
        }
            
        //* вычисления блоков по столбцу *//
        for (int b = 0; b < factor; ++b)
        {
            if (b == t)
                continue;

            // высота текущего (b-го по t-ому столбцу) блока
            const int ABheight = min(B, N - b*B);
            // указатель на начало текущего блока
            double* ABi = data + b*B*N + t*ABheight*B;
            for (int k = 0; k < CBsize; ++k)
            {
                const double* CBk = CB + k*CBsize;
                for (int i = 0; i < ABheight; ++i)
                {
                    // указатель на i-ую строку текущего блока
                    double* Ai = ABi + i * CBsize;
                    // к-ый элемент i-ой строки текущего блока
                    const double Aik = Ai[k];
                    for (int j = 0; j < CBsize; ++j)
                    {
                        Ai[j] = min(Ai[j], Aik + CBk[j]);
                    }
                }
            }
        }

        //* вычисления фазы 3 *//
        for (int b1 = 0; b1 < factor; ++b1)
        {
            if (b1 == t)
                continue;
                
            // высота текущей (b1-ой) полосы 
            const int ABheight = min(B, N - b1*B);
            // указатель на начало блока с координатами (b1, t)
            const double* b1_t_block = data + b1*B*N + t*ABheight*B;
            for (int b2 = 0; b2 < factor; ++b2)
            {
                if (b2 == t)
                    continue;

                // ширина текущего (b2-ого) столбца
                const int ABwidth = min(B, N - b2*B);
                // указатель на начало текущего блока (координаты = (b1,b2))
                double* ABi = data + b1*B*N + b2*ABheight*B;
                // указатель на начало блока с координатами (t,b2)
                const double* t_b2_block = data + t*B*N + b2*CBsize*B;
                for (int k = 0; k < CBsize; ++k)
                {
                    // указатель на k-ую строку блока с координатами (t, b2)
                    const double* t_b2_block_k_line = t_b2_block + k*ABwidth;
                    for (int i = 0; i < ABheight; ++i)
                    {
                        double* Ai = ABi + i * ABwidth;
                        const double Aik = b1_t_block[i*CBsize + k];
                        for (int j = 0; j < ABwidth; ++j)
                        {
                            Ai[j] = min(Ai[j], Aik + t_b2_block_k_line[j]);
                        }
                    }
                }
            }
        }
    }
    return data;
}
