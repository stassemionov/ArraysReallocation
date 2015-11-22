#include "service.h"

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>

using std::setw;

// M x N matrix generating
double* generate(double* A, const int M, const int N, const double lbound, const double ubound)
{
	double rnd = (ubound - lbound) / RAND_MAX;
	srand((int)time(NULL));
	for (long i = 0; i < M*N; ++i)
	{
		A[i] = lbound + rand() * rnd;
	}
	return A;
}

double* simple_fill(double* A, const int M, const int N)
{
	for (long i = 0; i < M*N; ++i)
	{
		A[i] = i+1;
	}
	return A;
}

void print_to(ostream& ostr, const double* A, const int M, const int N, const int place)
{
	const double* a = A;
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

double compare_arrays(const double* A, const double* B, const long long N)
{
	double s = 0;
	for (long i = 0; i < N; ++i)
	{
		s += abs(A[i] - B[i]);
	}
	return s;
}