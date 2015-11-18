#ifndef _SERVICE_H_
#define _SERVICE_H_

#include <iostream>

using namespace std;

// Generate M x N matrix. Elements are between lbound and ubound
double* generate(double* A, const int M, const int N, const double lbound = 0, const double ubound = 1);

double* simple_fill(double* A, const int M, const int N);

// Out standart-layout M x N matrix to the stream ostr
void print_to(ostream& ostr, const double* A, const int M, const int N, const int place = 10);
/*
// Out block-layout M x N matrix to the stream ostr with block highlighting
void block_print_to(ostream& ostr, double* A, int M, int N, int m, int n, int place = 10);
*/

double compare_arrays(const double* A, const double* B, const long long N);

#endif _SERVICE_H_