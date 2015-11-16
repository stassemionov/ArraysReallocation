#ifndef _SERVICE_H_
#define _SERVICE_H_

#include <iostream>

using namespace std;

// Generate M x N matrix. Elements are between lbound and ubound
double* generate(double* A, int M, int N, double lbound = 0, double ubound = 1);

double* simple_fill(double* A, int M, int N);

// Out standart-layout M x N matrix to the stream ostr
void print_to(ostream& ostr, double* A, int M, int N, int place = 10);

// Out block-layout M x N matrix to the stream ostr with block highlighting
void block_print_to(ostream& ostr, double* A, int M, int N, int m, int n, int place = 10);

double compare_arrays(double* A, double* B, long N);

#endif _SERVICE_H_