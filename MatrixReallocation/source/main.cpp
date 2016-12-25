#include "reallocation.h"
#include "service.h"
#include "multiplication.h"
#include "testing.h"

//#include <fstream>
//#include <ctime>

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

/*
C = A * B
Размещение матрицы A
- Матрица A разбивается на блоки, блоки - на подблоки
- блоки размещены по строкам
- подблоки размещены по столбцам
- элементы матрицы внутри подблока размещены по столбцам

Размещение матрицы B и C
- матрицы разбиваются на блоки
- блоки размещены по строкам.
- элементы внутри блоков размещены по строкам
*/


int main()
{
    setlocale(LC_ALL, "");

    InitDispatchSystem();

//    const TaskClass& floyd_params   = read_floyd_algorythm_parameters();
//    const TaskClass& qr_params      = read_qr_parameters();
//    const pair<TaskClass, TaskClass> mult_params = read_multiplication_parameters();
//    const TaskClass& realloc_params = read_reallocation_test_parameters();

//    floyd_test(floyd_params, true);
//    qralg_test(qr_params, true);
//    matrix_multiplication_tests(mult_params.first, mult_params.second, true);
//    reallocation_test(realloc_params, false);

    double* A = new double[25];
    for (int i = 0; i < 25; ++i)
    {
        A[i] = i;
    }
    standard_to_transposed_double_block_layout_reallocation(A, 5, 5, 5, 5, 3, 4);
    for (int i = 0; i < 25; ++i)
    {
        printf("%lf\n", A[i]);
    }
    printf("\n");
    delete[] A;

    TurnOffDispatchSystem();
}