#include "service.h"

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>

using namespace std;

#define _MIN_(x,y) (((x) < (y)) ? (x) : (y))

// M x N matrix generating
double* generate(double* A, int M, int N, double lbound, double ubound)
{
	double rnd = (ubound - lbound) / RAND_MAX;
	srand((int)time(NULL));
	for (long i = 0; i < M*N; ++i)
	{
		A[i] = lbound + rand() * rnd;
	}
	return A;
}

double* simple_fill(double* A, int M, int N)
{
	for (long i = 0; i < M*N; ++i)
	{
		A[i] = i+1;
	}
	return A;
}

void print_to(ostream& ostr, double* A, int M, int N, int place)
{
	double* a = A;
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			ostr << setw(place) << a[j];
		}
		a += N;
		ostr << '\n';
	}
	ostr << '\n';
}

// Block-by-block output of block-allocated array
void block_print_to(ostream& ostr, double* A, int M, int N, int m, int n, int place)
{/*
	for (int i = 0; i < M*N; )
	{
		for (int k = 0; k < _MIN_(m*n,M*N%m*n); ++j)
		{
			ostr << setw(place) << A[i];
		}
		ostr << '\n';
	}
	ostr << '\n';*/
}

double compare_arrays(double* A, double* B, long N)
{
	double s = 0;
	for (long i = 0; i < N; ++i)
	{
		s += abs(A[i] - B[i]);
	}
	return s;
}