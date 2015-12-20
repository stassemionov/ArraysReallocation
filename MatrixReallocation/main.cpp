#include <cstdlib>
#include <fstream>

#include "omp.h"

#include "reallocation.h"
#include "service.h"

using std::cout;
using std::ifstream;

int main()
{
    setlocale(LC_ALL, "");
    ostream& ostr = cout;
    ifstream file;
    file.open("../params.txt");
    // File needs to have a structure as follows (by a rows):
    // <number of rows    in matrix>\n          // N1
    // <number of columns in matrix>\n          // N2
    // <number of rows    in main block>\n      // b1
    // <number of columns in main block>\n      // b2
    // <number of rows    in small block>\n     // db1
    // <number of columns in small block>       // db2
    int N1, N2, b1, b2, db1, db2;
    file >> N1 >> N2 >> b1 >> b2 >> db1 >> db2;
    file.close();

    ostr << "N1  = " << N1 << "\n";
    ostr << "N2  = " << N2 << "\n";
    ostr << "b1  = " << b1 << "\n";
    ostr << "b2  = " << b2 << "\n";
    ostr << "db1 = " << db1 << "\n";
    ostr << "db2 = " << db2 << "\n\n";

    double* data = new double[N1*N2];
    ostr << "Заполнение массива...     ";
    double time_ = omp_get_wtime();
    simple_fill(data, N1, N2);
    time_ = omp_get_wtime() - time_;
    ostr << time_ << " секунд\n";
    print_to(ostr, data, N1, N2, 4);

    double* data_copy = new double[N1*N2];
    memcpy(data_copy, data, N1*N2*sizeof(double));

    ostr << "Вычисление эталона...     ";
    time_ = omp_get_wtime();
    double* dt = get_reallocated(data, N1, N2, b1, b2);
    time_ = omp_get_wtime() - time_;
    ostr << time_ << " секунд\n";
//  print_to(ostr, dt, N1, N2, 4);

    ostr << "Блочное переразмещение... ";
    time_ = omp_get_wtime();
    standard_to_block_layout_reallocation(data, N1, N2, b1, b2);
    time_ = omp_get_wtime() - time_;
    ostr << time_ << " секунд\n\n";
//  print_to(ostr, data, N1, N2, 4);

    ostr << "Норма разности: ";
    ostr << compare_arrays(dt, data, N1*N2) << "\n\n";

    ostr << "Двойное блочное переразмещение... ";
    time_ = omp_get_wtime();
    standard_to_double_block_layout_reallocation(data_copy, N1, N2, b1, b2, db1, db2);
    time_ = omp_get_wtime() - time_;
    ostr << time_ << " секунд\n\n";
    print_to(ostr, data_copy, N1, N2, 4);
}
