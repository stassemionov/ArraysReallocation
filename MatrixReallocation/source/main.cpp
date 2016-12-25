#include "reallocation.h"
#include "service.h"
#include "multiplication.h"
#include "testing.h"

//#include <fstream>
#include <ctime>

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
//    reallocation_test(realloc_params, false);

    int N1 = 4014, N2 = 4033, B1 = 411, B2 = 422, D1 = 86, D2 = 85;
    double* A = new double[N1 * N2], *B = new double[N1 * N2];
    for (int i = 0; i < N1 * N2; ++i)
    {
        A[i] = i;
    }
    memcpy(B, A, N1 * N2 * sizeof(double));

    double time_ = clock();

    InitDispatchSystem();
    standard_to_transposed_double_block_layout_reallocation(A, N1, N2, B1, B2, D1, D2);
    printf("\n AAAA \n");
    transposed_double_block_to_standard_layout_reallocation(A, N1, N2, B1, B2, D1, D2);
    printf("\n %lf %lf \n",
        compare_arrays(A, B, N1 * N2),
        (clock() - time_) / (1.0 * CLOCKS_PER_SEC));


    time_ = clock();

    standard_to_double_block_layout_reallocation(A, N1, N2, B1, B2, D1, D2);
    printf("\n AAAA \n");
    double_block_to_standard_layout_reallocation(A, N1, N2, B1, B2, D1, D2);
    TurnOffDispatchSystem();
    printf("\n %lf %lf \n",
        compare_arrays(A, B, N1 * N2),
        (clock() - time_) / (1.0 * CLOCKS_PER_SEC));

    delete[] A;
    delete[] B;
}