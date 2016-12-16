#include "reallocation.h"
#include "service.h"
#include "multiplication.h"
#include "testing.h"
#include "testing.h"

#include <fstream>
#include <ctime>

// Сделать!
//  1) Найти нестандартные эффективные способы обхода циклов                             (СЛОЖНО)
//  2) Переписать комментарии в соответствии с последними правками                       (ПРОСТО)
//  3) Сделать алгоритм двойного блочного переразмещения с использованием
//     соответствующей индексной функции, а не через композицию двух блочных
//     переразмещений, сравнить производительность и оставить лучший                     (СЛОЖНО)
//  4) В разделе тестирования одинаковый код вынести в отдельные функции                 (ПРОСТО)
//  5) Пересчитать сложность основного алгоритма переразмещения                          (СРЕДНЕ)
//  6) Применить исключения для обработки ошибок                                         (ПРОСТО)
//  7) Сделать функцию тестирующую всё                                                   (ПРОСТО)
//  8) Сравнивать блочные с блочными, двойные блочные с двойными блочными                (ПРОСТО)
//  9) Изменить проверку для теста переразмещения                                        (ПРОСТО)

int main()
{
    setlocale(LC_ALL, "");

    InitDispatchSystem();

    const TaskClass& floyd_params   = read_floyd_algorythm_parameters();
    const TaskClass& qr_params      = read_qr_parameters();
    const pair<TaskClass, TaskClass> mult_params = read_multiplication_parameters();
//    const TaskClass& realloc_params = read_reallocation_test_parameters();

    TurnOffDispatchSystem();

    floyd_test(floyd_params, true);
    qralg_test(qr_params, true);
    matrix_multiplication_tests(mult_params.first, mult_params.second, true);
//    reallocation_test(realloc_params, false);
}