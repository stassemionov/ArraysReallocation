#include "reallocation.h"
#include "service.h"

#include <fstream>

#include "omp.h"

using std::cout;
using std::ifstream;

int main()
{
    setlocale(LC_ALL, "");
    ostream& ostr = cout;
    ifstream file;
    file.open("../params.txt");
    // Sizes data file needs to have structure as follows:
    // <number of rows    in matrix>        ['\n' | ' ']+
    // <number of columns in matrix>        ['\n' | ' ']+
    // <number of rows    in main block>    ['\n' | ' ']+
    // <number of columns in main block>    ['\n' | ' ']+
    // <number of rows    in small block>   ['\n' | ' ']+
    // <number of columns in small block>   [' ']+ <anything>
    int N1, N2, b1, b2, db1, db2;
    file >> N1 >> N2 >> b1 >> b2 >> db1 >> db2;
    file.close();

    ostr << "N1  = " << N1 << "\n";
    ostr << "N2  = " << N2 << "\n";
    ostr << "b1  = " << b1 << "\n";
    ostr << "b2  = " << b2 << "\n";
    ostr << "db1 = " << db1 << "\n";
    ostr << "db2 = " << db2 << "\n\n";

    TaskClass parameters;
    parameters.makeData(N1, N2, b1, b2, db1, db2);

    double* data = new double[N1*N2];
    ostr << "Заполнение массива...     ";
    double time_ = omp_get_wtime();
    simple_fill(data, N1, N2);
    time_ = omp_get_wtime() - time_;
    ostr << time_ << " секунд\n";
//  print_to(ostr, data, N1, N2, 4);

    double* data_copy = new double[N1*N2];
    memcpy(data_copy, data, N1*N2*sizeof(double));

    ostr << "Вычисление эталона...     ";
    time_ = omp_get_wtime();
    double* dt = standard_to_block_layout_reallocation_buf(data, parameters);
    time_ = omp_get_wtime() - time_;
    ostr << time_ << " секунд\n";
//  print_to(ostr, dt, N1, N2, 4);

    ostr << "Блочное переразмещение... ";
    time_ = omp_get_wtime();
    standard_to_block_layout_reallocation(data, parameters);
    time_ = omp_get_wtime() - time_;
    ostr << time_ << " секунд\n\n";
//  print_to(ostr, data, N1, N2, 4);

    ostr << "Норма разности: ";
    ostr << compare_arrays(dt, data, N1*N2) << "\n\n";

    ostr << "Вычисление эталона...     ";
    time_ = omp_get_wtime();
    double* ddt = standard_to_double_block_layout_reallocation_buf(
        data_copy, parameters);
    time_ = omp_get_wtime() - time_;
    ostr << time_ << " секунд\n";
//  print_to(ostr, ddt, N1, N2, 4);

    ostr << "Двойное блочное переразмещение... ";
    time_ = omp_get_wtime();
    standard_to_double_block_layout_reallocation(data_copy, parameters);
    time_ = omp_get_wtime() - time_;
    ostr << time_ << " секунд\n\n";
//  print_to(ostr, data_copy, N1, N2, 4);

    ostr << "Норма разности: ";
    ostr << compare_arrays(ddt, data_copy, N1*N2) << "\n\n";
}
