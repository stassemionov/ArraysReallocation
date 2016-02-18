#include "reallocation.h"
#include "service.h"
#include "multiplication.h"
#include "testing.h"

#include <fstream>
#include <ctime>

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
    // <number of columns in matrix>        ['\n' | ' ']+//////////////////////////////////////
    // <number of rows    in main block>    ['\n' | ' ']+
    // <number of columns in main block>    ['\n' | ' ']+
    // <number of rows    in small block>   ['\n' | ' ']+
    // <number of columns in small block>   [' ']+ <anything>
    int N1, N2, N3, b1, b2, db1, db2;
    file >> N1 >> N2 >> N3 >> b1 >> b2 >> db1 >> db2;
    file.close();

    ostr << "\n [> Параметры задачи:\n";
    ostr << "  > N1  = " << N1 << "\n";
    ostr << "  > N2  = " << N2 << "\n";
    ostr << "  > N3  = " << N3 << "\n";
    ostr << "  > b1  = " << b1 << "\n";
    ostr << "  > b2  = " << b2 << "\n";
    ostr << "  > db1 = " << db1 << "\n";
    ostr << "  > db2 = " << db2 << "\n\n";

    ostr << " [> Подбор параметров двойного блочного размещения...\n";
    Blocks optimal = select_optimal_double_block_size(1000,
                                                      15, 150, 1,
                                                      10, 40, 1);
    ostr << "    Результат: " << optimal.main_block << " и " << optimal.small_block << std::endl;

/*    TaskClass params_left, params_right, params_gen;
    params_left.makeData(N1, N2, b1, b2, db1, db2);
    params_right.makeData(N2, N3, b2, b1, db2, db1);
    params_gen.makeData(N1, N3, b1, b1, db1, db1);

    double* left_mat = new double[N1*N2];
    double* rgt_mat  = new double[N2*N3];
    double* left_mat_copy = new double[N1*N2];
    double* rgt_mat_copy = new double[N2*N3];
    double* gen_mat1 = new double[N1*N3];
    double* gen_mat2 = new double[N1*N3];
    double* gen_mat3 = new double[N1*N3];
//    memset(gen_mat1, 0, N1*N3*sizeof(double));
//    memset(gen_mat2, 0, N1*N3*sizeof(double));

    const int memory_val = static_cast<int>( ceil((8.0 * (N1*N2 * 2 + N2*N3 * 2 + N1*N3 * 3)) / (1024 * 1024)) );
    ostr << "\n [> Объем выделенной памяти:                       ";
    ostr << "[ " << memory_val << " Мб ]\n";

    ostr << "\n [> Заполнение массивов...                         ";
    double time_ = clock();
    generate(left_mat, N1, N2);
    generate(rgt_mat,  N2, N3);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд ]\n";

    ostr << "\n [> Умножение матриц (стандартный алгоритм)...     ";
    time_ = clock();
    matrix_multiplication_standard(gen_mat1, left_mat, rgt_mat,
                                   N1, N2, N3);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд ]\n";

    ostr << "\n [> Умножение матриц (тайлинг)...                  ";
    time_ = clock();
    matrix_multiplication_tiled(gen_mat2, left_mat, rgt_mat,
                                params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд ]\n";

    ostr << "  > Норма разности: ";
    ostr << compare_arrays(gen_mat1, gen_mat2, N1*N3) << "\n";
    memset(gen_mat2, 0, N1*N3*sizeof(double));

    ostr << "\n [> Умножение матриц (двойной тайлинг)...          ";
    time_ = clock();
    matrix_multiplication_double_tiled(gen_mat2, left_mat, rgt_mat,
        params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд ]\n";

    ostr << "  > Норма разности: ";
    ostr << compare_arrays(gen_mat1, gen_mat2, N1*N3) << "\n";
    memset(gen_mat2, 0, N1*N3*sizeof(double));


    
    memcpy(left_mat_copy, left_mat, N1*N2*sizeof(double));
    ostr << "\n [> Блочное переразмещение...\n";
    ostr << "  > Левая матрица...                               ";
    time_ = clock();
    standard_to_block_layout_reallocation(left_mat, params_left);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд ]\n";

    memcpy(rgt_mat_copy, rgt_mat, N2*N3*sizeof(double));
    ostr << "  > Правая матрица...                              ";
    time_ = clock();
    standard_to_block_layout_reallocation(rgt_mat, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд ]\n";

    ostr << "\n [> Блочное умножение матриц (тайлинг)...          ";
    time_ = clock();
    block_matrix_multiplication_tiled(gen_mat2, left_mat, rgt_mat,
                                      params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд ]\n";
    
    memcpy(gen_mat3, gen_mat1, N1*N3*sizeof(double));
    ostr << "\n [> Вычисление эталона размещения...               ";
    time_ = clock();
    standard_to_block_layout_reallocation(gen_mat1, params_gen);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд ]\n";

    ostr << "  > Норма разности: ";
    ostr << compare_arrays(gen_mat1, gen_mat2, N1*N3) << "\n";
    memset(gen_mat2, 0, N1*N3*sizeof(double));
    



    ostr << "\n [> Двойное блочное переразмещение...\n";
    ostr << "  > Левая матрица...                               ";
    time_ = clock();
    standard_to_double_block_layout_reallocation(left_mat_copy, params_left);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд]\n";

    ostr << "  > Правая матрица...                              ";
    time_ = clock();
    standard_to_double_block_layout_reallocation(rgt_mat_copy, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд]\n";

    ostr << "\n [> Блочное умножение матриц (двойной тайлинг)...  ";
    time_ = clock();
    block_matrix_multiplication_double_tiled(gen_mat2, left_mat_copy, rgt_mat_copy,
        params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд]\n";

    ostr << "\n [> Вычисление эталона размещения...               ";
    time_ = clock();
    standard_to_double_block_layout_reallocation(gen_mat3, params_gen);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    ostr << "[ " << time_ << " секунд]\n";

    ostr << "  > Норма разности: ";
    ostr << compare_arrays(gen_mat3, gen_mat2, N1*N3) << "\n\n";
    */
}
