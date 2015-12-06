#ifndef _SERVICE_H_
#define _SERVICE_H_

#include <iostream>

using std::ostream;

// Generate M x N matrix. Elements are between lbound and ubound
double* generate(double* data_ptr, const int rows_count, const int cols_count,
    const double lbound = 0, const double ubound = 1);

// Matrix is filling with nubers [1, ..., M*N] with indices ascending order
double* simple_fill(double* data_ptr,
    const int rows_count, const int cols_count);

// Out standart-layout M x N matrix to the stream ostr
void print_to(ostream& ostr, const double* data_ptr, const int rows_count,
    const int cols_count, const int place = 10);

// Use sum of absolute values of arrays elements differences
// as measure of arrays difference
double compare_arrays(const double* data1, const double* data2, const int len);

#endif  // _SERVICE_H_
