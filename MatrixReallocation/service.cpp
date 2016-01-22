#include "service.h"

#include <iomanip>

using std::setw;
using std::vector;

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

bool m_find(const int& val, const vector<int>& vec)
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

    int l = 0;
    int r = static_cast<int>(vec.size()) - 1;
    int m = 0;
    while (l < r)
    {
        m = (l + r) / 2;
        if (vec[m] < val)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    return vec[r] == val;
}
