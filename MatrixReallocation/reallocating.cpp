#include "reallocating.h"

#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

int MM, NN, b1, b2;

inline long long f_ind(const long long& i);

double* get_reallocated(const double* A, const int M, const int N, const int m, const int n)
{
	MM = M;	// чтобы каждый раз не передавать эти параметры в индексную функцию
	NN = N;
	b1 = m;
	b2 = n;
	double* result = new double[M*N];
	for (long long i = 0; i < M*N; ++i)
	{
		result[f_ind(i)] = A[i];
	}
	return result;
}

inline long long f_ind(const long long& i)
{
	int&& row_i = (int) (i / NN);
	int&& col_i = i % NN;
	int&& row_b = row_i / b1;	// координаты блока
	int&& col_b = col_i / b2;

	int&& M_b = (int) ceil((double) MM / b1);	// число блоков по столбцу (в направлении столбца)
	int&& N_b = (int) ceil((double) NN / b2);	// число блоков по строке (в направлении строки)
	int&& I = row_b * N_b + col_b;	// номер блока

//	int row_i_loc = row_i % b1;	// координаты элемента в блоке
//	int col_i_loc = col_i % b2;

	int w = ( (col_b == N_b - 1) && (NN % b2 != 0) ) ? NN % b2 : b2;
	int&& shift = (row_i % b1) * w + (col_i % b2);
	int&& dif = NN % b2;
	if ((row_b == M_b - 1) && (MM % b1 != 0))	// для нижних блоков нужно вычислить смещение границы сетки от границы матрицы
	{		
		if (dif == 0)
			return	I*b1*b2 + shift - col_b*b2*(b1 - MM % b1);
		else
			return	I*b1*b2 + shift - row_b*b1*(b2 - dif) - col_b*b2*(b1 - MM % b1);
	}
	else
	{
		if (dif == 0)
			return I*b1*b2 + shift;
		else
			return	I*b1*b2 + shift - row_b*b1*(b2 - dif);
	}
}

bool is_new_cycle(const long long i, const vector<long long>& sdr_vec)
{
	// подсчитывать длину цикла? или другие его параметры?
	long long next = i;
	do
	{
		for (size_t j = 0; j < sdr_vec.size(); ++j)
		{
			if (next == sdr_vec[j])
			{
				return false;
			}
		}
		next = f_ind(next);
	} 
	while (next != i);
	return true;
}

double* case1(double* A, const int M, const int N, const int m, const int n)
{
	vector<long long> sdr_vec; // system of distinct representatives of cycles

	// * Обучение алгоритма текущей картине распределения циклов
	long long i = n;	// Итерации начинаются с элемента 'n', т.к. первые 'n' элементов не требуют перемещения.
	long long max_index = i;	// максимальный индекс среди элементов, попавших в очередной цикл
	const long long iteration_count = m * N - 2 * n;
	long long it = 0;

	while (it < iteration_count)
	{
		const long long first_in_cycle = i;
		sdr_vec.push_back(first_in_cycle);
		do
		{
			i = f_ind(i);
			if (max_index < i)
			{
				max_index = i;
			}
			it += n;
		} 
		while (first_in_cycle != i);	// пока не попадем в начало цикла (элемент, с которого начинали)

		if (it < iteration_count)
		{
			long long next_cycle_begining = max_index - n;
			while (!is_new_cycle(next_cycle_begining, sdr_vec))	// поиск следующего цикла
			{
				next_cycle_begining -= n;
			}
			i = next_cycle_begining;
			max_index = i;
		}		
	}
	// * Конец обучения
		
	double* s = new double[n];
	double* b = new double[n];
	const int len = n*sizeof(double);
	const int stripe_size = m*N;

	for (int v = 0; v < M / m; ++v)	// по горизонтальным полосам
	{
		size_t cycle_counter = 0;		
		while (cycle_counter < sdr_vec.size())	// действия в пределах выбранной полосы
		{
			i = sdr_vec[cycle_counter++] + v*stripe_size;
			memcpy(s, A + i, len);
			const long long first_in_cycle = i;
			do
			{
				i = f_ind(i);
			
				memcpy(b, A + i, len);
				memcpy(A + i, s, len);
				memcpy(s, b, len);
			} 
			while (first_in_cycle != i);	// пока не попадем в начало цикла (элемент, с которого начинали)
		}
	}
/*
	const int&& last_stride_h = M % b1;
	if (last_stride_h > 1)
	{
		long i = (M - last_stride_h) * N + n;
		memcpy(s, A + i, len);
		for (long k = 0; k < last_stride_h*N - 2*n; k += n)
		{
			i = f_ind(i);

			memcpy(b, A + i, len);
			memcpy(A + i, s, len);
			memcpy(s, b, len);
		}
	}*/

	delete[] s;
	delete[] b; 

	return A;
}

double* block_reallocate_matrix(double* A, const int M, const int N, const int m, const int n)
{
	MM = M;	// чтобы каждый раз не передавать эти параметры в индексную функцию
	NN = N;
	b1 = m;
	b2 = n;
	if (N % n == 0)
	{
		return case1(A, M, N, m, n);
	}
/*	else
	{
		;
		return case1(A, M, N, m, n);
	}*/
	return A;
}
