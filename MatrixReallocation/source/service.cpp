#include "service.h"

#include <iomanip>
#include <fstream>
#include <cstring>
#include <algorithm>

using std::min;
using std::max;
using std::setw;
using std::ifstream;
using std::make_pair;

double* generate(double* data_ptr, 
    const int rows_count, const int cols_count,
    const double lbound, const double ubound)
{
    double rnd = (ubound - lbound) / RAND_MAX;
    srand(static_cast<int>(time(NULL)));
    const size_t length = static_cast<size_t>(rows_count * cols_count);
    for (size_t i = 0; i < length; ++i)
    {
        data_ptr[i] = lbound + rand() * rnd;
    }
    return data_ptr;
}

double* simple_fill(double* data_ptr,
    const int rows_count, const int cols_count)
{
    const size_t length = static_cast<size_t>(rows_count * cols_count);
    for(size_t i = 0; i < length; ++i)
    {
        data_ptr[i] = (double) i;
    }
    return data_ptr;
}

void print_to(ostream& ostr, const double* data_ptr,
    const int rows_count, const int cols_count, const int place)
{
    const double* row = data_ptr;
    for (int i = 0; i < rows_count; ++i)
    {
        for (int j = 0; j < cols_count; ++j)
        {
            ostr << setw(place) << row[j];
        }
        row += cols_count;
        ostr << '\n';
    }
    ostr << '\n';
}

double compare_arrays(
    const double* data1, const double* data2, const size_t length)
{
    double s(0.0);
    for (size_t i = 0; i < length; ++i)
    {
        s += fabs(data1[i] - data2[i]);
    }
    return s;
}

int gcd(const int u, const int v)
{
    // simple cases (termination)
    if (u == v)
        return u;

    if (u == 0)
        return v;

    if (v == 0)
        return u;

    // look for factors of 2
    if (~u & 1)     // u is even
    {
        if (v & 1)  // v is odd
            return gcd(u >> 1, v);
        else    // both u and v are even
            return gcd(u >> 1, v >> 1) << 1;
    }

    if (~v & 1)     // u is odd, v is even
        return gcd(u, v >> 1);

    // reduce larger argument
    if (u > v)
        return gcd((u - v) >> 1, v);

    return gcd((v - u) >> 1, u);
}

int bin_search(const int val, const vector<int>& vec)
{
    if (vec.empty())
    {
        return false;
    }
    if (val > vec.back())
    {
        return false;
    }
    if (val < vec.front())
    {
        return false;
    }

    int l(0);
    int r = static_cast<int>(vec.size()) - 1;
    int m(0);
    while (l < r)
    {
        m = (l + r) >> 1;
        if (vec[m] < val)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    return (vec[r] == val) ? r : -1;
}

pair<TaskClass, TaskClass> read_multiplication_parameters(
    const string& file_name)
{
    // Task paramers file needs to have structure as follows:
    // <number of rows         in left       matrix>                  ['\n' | ' ']+
    // <number of columns/rows in left/rigth matrix>                  ['\n' | ' ']+
    // <number of columns      in right      matrix>                  ['\n' | ' ']+
    // <number of rows         in main block of left       matrix>    ['\n' | ' ']+
    // <number of columns/rows in main block of left/right matrix>    ['\n' | ' ']+
    // <number of columns      in main block of right      matrix>    ['\n' | ' ']+
    // <number of rows         in small block of left       matrix>   ['\n' | ' ']+
    // <number of columns/rows in small block of left/right matrix>   ['\n' | ' ']+
    // <number of rows         in small block of right      matrix>   ['\n' | ' ']+ <anything>

    ifstream file;
    file.open(file_name);

    int N1, N2, N3, b1, b2, b3, db1, db2, db3;
    file >> N1 >> N2 >> N3 >> b1 >> b2 >> b3 >> db1 >> db2 >> db3;
    file.close();

    return make_pair<TaskClass, TaskClass>(
        TaskClass(N1, N2, b1, b2, db1, db2),
        TaskClass(N2, N3, b2, b3, db2, db3));
}

TaskClass read_reallocation_test_parameters(const string& file_name)
{
    // Task paramers file needs to have structure as follows:
    // <number of rows    in matrix>                  ['\n' | ' ']+
    // <number of columns in matrix>                  ['\n' | ' ']+
    // <number of rows    in main block of matrix>    ['\n' | ' ']+
    // <number of columns in main block of matrix>    ['\n' | ' ']+
    // <number of rows    in small block of matrix>   ['\n' | ' ']+
    // <number of rows    in small block of matrix>   ['\n' | ' ']+ <anything>

    ifstream file;
    file.open(file_name);

    int N1, N2, b1, b2, db1, db2;
    file >> N1 >> N2 >> b1 >> b2 >> db1 >> db2;
    file.close();

    return TaskClass(N1, N2, b1, b2, db1, db2);
}

TaskClass read_floyd_algorythm_parameters(const string& file_name)
{
    // Task paramers file needs to have structure as follows:
    // <number of rows/columns in matrix>   ['\n' | ' ']+
    // <number of rows in block of matrix>  ['\n' | ' ']+ <anything>
    
    ifstream file;
    file.open(file_name);

    int N, b;
    file >> N >> b;
    file.close();

    return TaskClass(N, N, b, b);
}

TaskClass read_qr_parameters(const string& file_name)
{
    // Task paramers file needs to have structure as follows:
    // <number of rows/columns in matrix>                  ['\n' | ' ']+
    // <number of rows         in main block of matrix>    ['\n' | ' ']+
    // <number of columns      in main block of matrix>    ['\n' | ' ']+
    // <number of rows         in small block of matrix>   ['\n' | ' ']+
    // <number of columns      in small block of matrix>   ['\n' | ' ']+ <anything>

    ifstream file;
    file.open(file_name);

    int N, b1, b2, db1, db2;
    file >> N >> b1 >> b2 >> db1 >> db2;
    file.close();

    return TaskClass(N, N, b1, b2, db1, db2);
}
