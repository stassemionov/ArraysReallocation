﻿#ifndef _SERVICE_H_
#define _SERVICE_H_

#include "taskdata.h"

#include <iostream>
#include <vector>

using std::ostream;
using std::vector;
using std::string;
using std::pair;

// Generate M x N matrix. Elements are between lbound and ubound
double* generate(double* data_ptr, const int rows_count, const int cols_count,
    const double lbound = 0, const double ubound = 1);

// Matrix is filling with nubers [1, ..., M*N] with indices ascending order
double* simple_fill(double* data_ptr,
    const int rows_count, const int cols_count);

// Matrix data output to 'ostr' stream
void print_to(ostream& ostr, const double* data_ptr, const int rows_count,
    const int cols_count, const int place = 10);

// Use sum of absolute values of arrays elements differences
// as measure of arrays difference
double compare_arrays(const double* data1, const double* data2,
                      const size_t length);

// Binary Euclidean agorythm
int gcd(const int u, const int v);

// Search 'val' in vector 'vec' with binary search,
// return true or false
int bin_search(const int val, const vector<int>& vec);

// Writes task parameters from file with name 'file_name'.
// Returns 3-element dynamic array of TaskClass objects,
// containing elements in following order:
// left matrix parameters, right matrix parameters, generator matrix parameters
pair<TaskClass, TaskClass> read_multiplication_parameters(
    const string& file_name = "../resource/parameter/mult_parameters.txt");

// Writes task parameters from file with name 'file_name'.
// Returns pointer to dynamic-allocated TaskClass object
TaskClass read_reallocation_test_parameters(
    const string& file_name = "../resource/parameter/realloc_parameters.txt");

// Writes task parameters from file with name 'file_name'.
// Returns pointer to dynamic-allocated TaskClass object
TaskClass read_floyd_algorythm_parameters(
    const string& file_name = "../resource/parameter/floyd_parameters.txt");

// Writes task parameters from file with name 'file_name'.
// Returns pointer to dynamic-allocated TaskClass object
TaskClass read_qr_parameters(
    const string& file_name = "../resource/parameter/qr_parameters.txt");

#endif  // _SERVICE_H_
