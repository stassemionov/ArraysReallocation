#include "reallocating.h"

#include <cstdlib>
#include <cmath>

using namespace std;

#include <iostream>
#include <iomanip>

int MM, NN, b1, b2;

inline long f_ind(const long& i)
{
	int row_i = i / NN;
	int col_i = i % NN;
	int row_b = row_i / b1;	// координаты блока
	int col_b = col_i / b2;

	int M_b = ceil((double) MM / b1);	// число блоков по столбцу (в направлении столбца)
	int N_b = ceil((double) NN / b2);	// число блоков по строке (в направлении строки)
	int I = row_b * N_b + col_b;	// номер блока

	int row_i_loc = row_i % b1;	// координаты элемента в блоке
	int col_i_loc = col_i % b2;

	int w = ( (col_b == N_b - 1) && (NN % b2 != 0) ) ? NN % b2 : b2;
	int shift = row_i_loc * w + col_i_loc;
	if (row_b == M_b - 1)	// для нижних блоков нужно вычислить смещение границы сетки от границей матрицы
	{
		if (MM % b1 != 0)
		{
			int&& dif = NN % b2;
			if (dif == 0)
				return	I*b1*b2 + shift - col_b*b2*(b1 - MM % b1);
			else
				return	I*b1*b2 + shift - row_b*b1*(b2 - dif) - col_b*b2*(b1 - MM % b1);
		}
		else
		{
			int&& dif = NN % b2;
			if (dif == 0)
				return	I*b1*b2 + shift;
			else
				return	I*b1*b2 + shift - row_b*b1*(b2 - NN % b2);
		}
	}
	else
	{
		int&& dif = NN % b2;
		if (dif == 0)
			return I*b1*b2 + shift;
		else
			return	I*b1*b2 + shift - row_b*b1*(b2 - dif);
	}
}

double* case1(double* A, int M, int N, int m, int n)
{
	long fi = 0;
	double* s = new double[b2];
	double* b = new double[b2];
	int len = b2*sizeof(double);

	for (int v = 0; v < M / b1; ++v)	// по горизонтальным полосам
	{
		long i = v * N + b2;	// Итерации начинаются с элемента b2, т.к. первые b2 элементов не требуют перемещения.
		memcpy(s, A + v*N + b2, len);
		for (long k = 0; k < b1*N - 2*b2; k += b2)	// действия в пределах выбранной полосы
		{
			fi = f_ind(i);

			memcpy(b, A + fi, len);
			memcpy(A + fi, s, len);
			memcpy(s, b, len);
		}
	}

	int L = ((M % b1 == 1) ? N + 2 * b2 : b2 + 1);

	return A;
}

double* block_reallocate_matrix(double* A, int M, int N, int m, int n)
{
	MM = M;
	NN = N;
	b1 = m;
	b2 = n;
	if (N % n == 0)
	{
		return case1(A, M, N, m, n);
	}
	else
	{
		return case1(A, M, N, m, n);
	}
}

double* buffer_reallocating(double* A, int M, int N, int m, int n)
{
	MM = M;
	NN = N;
	b1 = m;
	b2 = n;
	double* result = new double[M*N];
	for (long i = 0; i < M*N; ++i)
	{
		result[f_ind(i)] = A[i];
	}
	return result;
}

/*
double* block_reallocate_matrix(double* A, int M, int N, int m, int n)
{
	long i = b2;	// Итерации начинаются с элемента b2, т.к. первые b2 элементов не требуют перемещения.
	long fi = 0;
	double* s = new double[b2];
	double* b = new double[b2];
	int L = ((M % b1 == 1) ? N+2*b2 : b2+1);
	int len = b2*sizeof(double);

	memcpy(s, A + b2, len);

	for (long k = 0; k < M*N-L; k += b2)
	{
		fi = f_ind(i);
		//		cout << i + 1 << setw(3) << "-" << setw(3) << fi + 1 << endl;
		//		(i == fi) ? ++i : i = fi;

		if (i == fi)
		{
			i += b2;
			continue;
		}
		else
		{
			i = fi;
		}

		memcpy(b, A + fi, len);
		memcpy(A + fi, s, len);
		memcpy(s, b, len);

		//b = A[fi];
		//A[fi] = s;
		//s = b;
	}
	return A;
}
*/
