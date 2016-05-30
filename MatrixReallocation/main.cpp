#include "reallocation.h"
#include "service.h"
#include "multiplication.h"
#include "testing.h"

#include <fstream>
#include <ctime>

#include "qralg.h"

using std::ifstream;

// Сделать!
//  1) Сократить просмотр больших промежутков при поиске нового цикла                    (СЛОЖНО)
//  2) Определить более надежное условие параллельности циклов (да и это не подводит..)  (СРЕДНЕ)
//  3) Найти способ не заносить элементы параллельных циклов в help_vec                  (СРЕДНЕ)
//  4) Подумать и поэкспериментировать, можно ли,
//     имея sdr для [b1,b2,d1,d2] и для [b2,b3,d2,d3], получить sdr для [b1,b3,d1,d3] ?  (СЛОЖНО)
//  5) Применить экономичные структуры данных для оптимизации объема памяти
//     для хранения посещенных элементов                                                 (ПРОСТО)
//  6) Исправить ошибки в release-версиях функций переразмещения
//     (инициализация некоторых полей структур нулями)                                   (ПРОСТО)
//  7) Проход по циклу сделать с двух сторон (?)                                         (ПРОСТО)
//  8) Добавить тест для копирования куска матрицы                                       (ПРОСТО)
//  9) Изменить способ определения нового цикла                                          (ПРОСТО)
// 10) Разобраться с обратным переразмещением матрицы-генератора при умножении матриц    (СРЕДНЕ)
// 11) Сделать QR-разложение с двойным тайлингом                                         (СРЕДНЕ)

int main()
{
    setlocale(LC_ALL, "");

//    const TaskClass& floyd_params   = *read_floyd_algorythm_parameters();
    const TaskClass& qr_params      = *read_qr_parameters();
//    const TaskClass& mult_params    = *read_multiplication_parameters();
//    const TaskClass& realloc_params = *read_reallocation_test_parameters();
    
    qralg_test(qr_params, true);
}


/*
const TaskClass& params = *read_qr_parameters();
const int N = params.getDataRef().M_ROWS;

double* mat_t = new double[N*N];
double* mat_b = new double[N*N];
double* mat_db = new double[N*N];
double* mat_st = new double[N*N];

generate(mat_st, N, N);
memcpy(mat_b, mat_st, N*N*sizeof(double));
memcpy(mat_db, mat_st, N*N*sizeof(double));

int rowc = 1;
int colc = 23;

double time_ = clock();
// Обычный алгоритм
for (int i = rowc; i < N; ++i)
{
memset(mat_db + i*N + colc, 0,
(N - colc)*sizeof(double));
}
time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
std::cout << "STANDARD COMPLETE! " << time_ << std::endl << std::endl;

// Алгоритм для блочного размещения
const DoubleBlockReallocationInfo* info =
standard_to_double_block_layout_reallocation(mat_b, params);
time_ = clock();
filmat_part(mat_b, 0, params.getDataRef(), rowc, colc);
time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
std::cout << "FILMAT COMPLETE! " << time_ << std::endl << std::endl;
double_block_to_standard_layout_reallocation(mat_b, *info);

std::cout << "\n" << compare_arrays(mat_b, mat_db, N*N) << "\n";

//    print_to(std::cout, mat_db, N, N, 10);
//    print_to(std::cout, mat_b, N, N, 10);
*/