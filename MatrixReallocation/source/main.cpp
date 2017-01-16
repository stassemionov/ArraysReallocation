#include "reallocation.h"
#include "service.h"
#include "multiplication.h"
#include "testing.h"
#include "mappingfuncs.h"

#include <iostream>
#include <ctime>
#include <omp.h>

using std::cout;
using std::endl;

// Сделать!
//  1) Найти нестандартные эффективные способы обхода циклов                             (СЛОЖНО)
//  2) Переписать комментарии в соответствии с последними правками                       (ПРОСТО)
//  3) Сделать алгоритм двойного блочного переразмещения с использованием
//     соответствующей индексной функции, а не через композицию двух блочных
//     переразмещений, сравнить производительность и оставить лучший                     (СЛОЖНО)
//  4) В разделе тестирования одинаковый код вынести в отдельные функции                 (ПРОСТО)
//  5) Применить исключения для обработки ошибок                                         (ПРОСТО)
//  6) Сделать функцию тестирующую всё                                                   (ПРОСТО)
//  7) Изменить проверку для теста переразмещения                                        (ПРОСТО)

int main()
{
    setlocale(LC_ALL, "");

//    const TaskClass& floyd_params   = read_floyd_algorythm_parameters();
//    const TaskClass& qr_params      = read_qr_parameters();
//    const pair<TaskClass, TaskClass> mult_params = read_multiplication_parameters();
//    const TaskClass& realloc_params = read_reallocation_test_parameters();

//    floyd_test(floyd_params, true);
//    qralg_test(qr_params, true);
//    matrix_multiplication_tests(mult_params.first, mult_params.second, true);
//    reallocation_test(realloc_params, true);

    int N1 = 8192, N2 = N1,
        B1 = 512, B2 = 256,
        D1 = 8, D2 = 256;
    double
        *A = new double[N1 * N2],
        *B = new double[N1 * N2];
    for (int i = 0; i < N1 * N2; ++i)
    {
        A[i] = i + 1;
    }
//    memcpy(B, A, N1 * N2 * sizeof(double));
    map_with_transposed_double_block_layout(B, A, N1, N2, B1, B2, D1, D2);
//    print_to(cout, B, 1, N1*N2, 3);

//    double time_ = omp_get_wtime();
//
//    InitDispatchSystem();
//    standard_to_transposed_double_block_layout_reallocation(A, N1, N1, 512, 256, 4, 256);
//    printf("\n %lf \n", omp_get_wtime() - time_);
//    standard_to_double_block_layout_reallocation(A, N1, N1, 256, 8, 256, 8);
//    standard_to_double_block_layout_reallocation(A, N1, N1, 4, 8, 4, 8);
//    printf("\n %lf \n", omp_get_wtime() - time_);
//    transposed_double_block_to_standard_layout_reallocation(A, N1, N1, 512, 256, 4, 256);
//    printf("\n %lf \n", omp_get_wtime() - time_);
//    double_block_to_standard_layout_reallocation(A, N1, N1, 256, 8, 256, 8);
//    double_block_to_standard_layout_reallocation(A, N1, N1, 4, 8, 4, 8);
//    TurnOffDispatchSystem();
//    printf("\n %lf \n", omp_get_wtime() - time_);
    

    double time_ = omp_get_wtime();

    InitDispatchSystem();
    standard_to_transposed_double_block_layout_reallocation(A, N1, N2, B1, B2, D1, D2);
    time_ = omp_get_wtime() - time_;
    cout << compare_arrays(A, B, N1 * N2) << " " << time_ << endl;

//    time_ = omp_get_wtime();
//    transposed_double_block_to_standard_layout_reallocation(A, N1, N2, B1, B2, D1, D2);
//    time_ = omp_get_wtime() - time_;
//    printf("%lf %lf \n",
//        compare_arrays(A, B, N1 * N2),
//        time_);
    
//    time_ = omp_get_wtime();
//
//    standard_to_double_block_layout_reallocation(A, N1, N2, B1, B2, D1, D2);
//    double_block_to_standard_layout_reallocation(A, N1, N2, B1, B2, D1, D2);
//    TurnOffDispatchSystem();
//    time_ = omp_get_wtime() - time_;
//    printf("%lf %lf \n",
//        compare_arrays(A, B, N1 * N2),
//        time_);

    delete[] A;
    delete[] B;
}