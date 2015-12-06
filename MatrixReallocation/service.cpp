#include "service.h"

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>

using std::setw;

double* generate(double* data_ptr, const int rows_count, const int cols_count,
    const double lbound, const double ubound)
{
    double rnd = (ubound - lbound) / RAND_MAX;
    srand(static_cast<int>(time(NULL)));
    for (int i = 0; i < rows_count*cols_count; ++i)
    {
        data_ptr[i] = lbound + rand() * rnd;
    }
    return data_ptr;
}

double* simple_fill(double* data_ptr,
    const int rows_count, const int cols_count)
{
    for (int i = 0; i < rows_count*cols_count; ++i)
    {
        data_ptr[i] = i+1;
    }
    return data_ptr;
}

// Matrix data output to 'ostr' stream
void print_to(ostream& ostr, const double* data_ptr,
    const int rows_count, const int cols_count, const int place)
{
    const double* a = data_ptr;
    for (int i = 0; i < rows_count; ++i)
    {
        for (int j = 0; j < cols_count; ++j)
        {
            ostr << setw(place) << a[j];
        }
        a += cols_count;
        ostr << '\n';
    }
    ostr << '\n';
}

double compare_arrays(const double* data1, const double* data2, const int len)
{
    double s = 0;
    for (int i = 0; i < len; ++i)
    {
        s += abs(data1[i] - data2[i]);
    }
    return s;
}
