#include "../include/reallocation.h"
#include "../include/service.h"
#include "../include/multiplication.h"
#include "../include/testing.h"

#include <fstream>
#include <ctime>

// Сделать!
//  1) Сократить просмотр больших промежутков при поиске нового цикла                    (СЛОЖНО)
//  2) Найти нестандартные эффективные способы обхода циклов                             (СЛОЖНО)
//  3) Выработать новый подход к применению доп. векторов                                (СЛОЖНО)
//  4) Переписать комментарии в соответствии с последними правками                       (ПРОСТО)
//  5) Сделать алгоритм двойного блочного переразмещения с использованием
//     соответствующей индексной функции, а не через композицию двух блочных
//     переразмещений, сравнить производительность и оставить лучший                     (СЛОЖНО)
//  6) В разделе тестирования одинаковый код вынести в отдельные функции                 (ПРОСТО)
//  7) Пересчитать сложность основного алгоритма переразмещения                          (СРЕДНЕ)
//  8) Применить исключения для обработки ошибок                                         (ПРОСТО)
//  9) Сделать функцию тестирующую всё                                                   (ПРОСТО)
// 10) Сравнивать блочные с блочными, двойные блочные с двойными блочными                (ПРОСТО)

int main()
{
    setlocale(LC_ALL, "");

    InitDispatchSystem();

//    const TaskClass& floyd_params   = read_floyd_algorythm_parameters();
//    const TaskClass& qr_params      = read_qr_parameters();
//    const pair<TaskClass, TaskClass> mult_params = read_multiplication_parameters();
    const TaskClass& realloc_params = read_reallocation_test_parameters();

//    floyd_test(floyd_params, true);
//    qralg_test(qr_params, true);
//    matrix_multiplication_tests(mult_params.first, mult_params.second, true);
    reallocation_test(realloc_params, true);
}