#include "reallocation.h"
#include "service.h"
#include "multiplication.h"
#include "testing.h"

#include <fstream>
#include <ctime>

#include "qralg.h"

using std::ifstream;

// Сделать!
// 1) Сократить просмотр больших промежутков при поиске нового цикла                    (СЛОЖНО)
// 2) Определить более надежное условие параллельности циклов (да и это не подводит..)  (СРЕДНЕ)
// 3) Найти способ не заносить элементы параллельных циклов в help_vec (СРЕДНЕ)
// 6) Подумать и поэкспериментировать, можно ли,
//    имея sdr для [b1,b2,d1,d2] и для [b2,b3,d2,d3], получить sdr для [b1,b3,d1,d3] ?  (СЛОЖНО)
// 7) Применить экономичные структуры данных для оптимизации объема памяти
//    для хранения посещенных элементов                                                 (ПРОСТО)
// 8) Исправить ошибку в функции стандартного QR-разложения (сравнение с maple)         (СРЕДНЕ)
// 9) Исправить ошибки в release-версиях функций переразмещения
//    (инициализация некоторых полей структур нулями)                                   (ПРОСТО)

int main()
{
    setlocale(LC_ALL, "");

//    const TaskClass* floyd_params   = read_floyd_algorythm_parameters();
    const TaskClass* qr_params      = read_qr_parameters();
//    const TaskClass* mult_params    = read_multiplication_parameters();
//    const TaskClass* realloc_params = read_reallocation_test_parameters();

    /*int N = 5;
    int B1 = 5;
    int B2 = 2;
    int D1 = 2;
    int D2 = 2;*/
//    TaskClass params(N,N,B1,B2,D1,D2);

    const TaskClass& params = *qr_params;
    const int NN = params.getDataRef().M_ROWS * params.getDataRef().M_ROWS;

    double* mat_t = new double[NN];
    double* mat_b  = new double[NN];
    double* mat_db = new double[NN];
    double* mat_st = new double[NN];
    simple_fill(mat_st, params.getDataRef().M_ROWS, params.getDataRef().M_ROWS);
    mat_st[0] *= -1;
    mat_st[1] *= -1;
    mat_st[2] *= -1;
    mat_st[11] *= -1;
    mat_st[12] *= -1;
    mat_st[13] *= -1;
    mat_st[22] *= -1;
    mat_st[23] *= -1;
    mat_st[24] *= -1;
    memcpy(mat_t,  mat_st, NN*sizeof(double));
    memcpy(mat_b,  mat_st, NN*sizeof(double));
    memcpy(mat_db, mat_st, NN*sizeof(double));

    std::cout << std::endl;
    print_to(std::cout, mat_st, params.getDataRef().M_ROWS, params.getDataRef().M_ROWS);
    
    std::cout << "\nSTART!\n\n";
    
    double time_ = clock();
    QR_WY_standard(mat_st, params.getDataRef().M_ROWS, params.getDataRef().B_COLS);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    std::cout << "STANDARD COMPLETE! " << time_ << std::endl << std::endl;

    time_ = clock();
    QR_WY_tiled(mat_t, params.getDataRef());
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    std::cout << "TILED COMPLETE! " << time_ << std::endl << std::endl;
    
    const BlockReallocationInfo* info =
        standard_to_block_layout_reallocation(mat_b, params);
    std::cout << "REALLOC COMPLETE!\n";
    time_ = clock();
    QR_WY_block(mat_b, params.getDataRef());
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    std::cout << "BLOCK COMPLETE! " << time_ << std::endl;
    block_to_standard_layout_reallocation(mat_b, *info);
    std::cout << "INVERSE REALLOC COMPLETE!\n\n";
    
    const DoubleBlockReallocationInfo* info2 =
        standard_to_double_block_layout_reallocation(mat_db, params);
    std::cout << "REALLOC COMPLETE!\n";
    time_ = clock();
    QR_WY_double_block(mat_db, params.getDataRef());
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    std::cout << "DOUBLE BLOCK COMPLETE! " << time_ << std::endl;
    double_block_to_standard_layout_reallocation(mat_db, *info2);
    std::cout << "INVERSE REALLOC COMPLETE!\n\n";

    std::cout << " \n DIFFERENCE DOUBLE - BLOCK [ " << compare_arrays(mat_db, mat_b, NN)  << " ]";
    std::cout << " \n DIFFERENCE DOUBLE - TILED [ " << compare_arrays(mat_db, mat_t, NN)  << " ]";
    std::cout << " \n DIFFERENCE DOUBLE - STAND [ " << compare_arrays(mat_db, mat_st, NN) << " ]";
    std::cout << " \n DIFFERENCE BLOCK  - TILED [ " << compare_arrays(mat_b, mat_t, NN)   << " ]";
    std::cout << " \n DIFFERENCE BLOCK  - STAND [ " << compare_arrays(mat_b, mat_st, NN)  << " ]";
    std::cout << " \n DIFFERENCE TILED  - STAND [ " << compare_arrays(mat_t, mat_st, NN)  << " ]\n ";

    std::cout << std::endl;
    print_to(std::cout, mat_b, params.getDataRef().M_ROWS, params.getDataRef().M_ROWS, 12);

    delete[] mat_t;
    delete[] mat_b;
    delete[] mat_db;
    delete[] mat_st;
}
