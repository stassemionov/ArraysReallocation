#include "testing.h"
#include "reallocation.h"
#include "multiplication.h"
#include "floydalg.h"
#include "qralg.h"
#include "service.h"

#include <cstring>
#include <ctime>
#include <algorithm>

using std::min;

Blocks select_optimal_double_block_size_multiplication(const int N,
    const int start_block_val,
    const int end_block_val,
    const int selection_block_step,
    const int start_double_block_val,
    const int end_double_block_val,
    const int selection_double_block_step,
    const bool         console_info_output)
{
    FILE* log_file;
    fopen_s(&log_file, "../resource/log/double_block_selection_test_log.txt", "w+");

    fprintf(log_file, "\n [> Подбор параметров двойного блочного размещения...\n");
    fprintf(log_file, "  > Обозначения:\n");
    fprintf(log_file, "     N - размер матрицы\n");
    fprintf(log_file, "     B и D - размеры большого и малого блоков\n");
    fprintf(log_file, "     MT - продолжительность умножения матриц в секундах\n");
    fprintf(log_file, "     RT - продолжительность переразмещения матриц в секундах\n");
    fprintf(log_file, "     WT = MT + RT\n");
    fprintf(log_file, "  > Параметры подбора:\n     N = %d, B = [%d:%d:%d], D = [%d:%d:%d]\n",
            N, start_block_val, end_block_val, selection_block_step,
            start_double_block_val, end_double_block_val, selection_double_block_step);

    if (console_info_output)
    {
        printf("\n [> Подбор параметров двойного блочного размещения...\n");
        printf("  > N = %d, B = [%d:%d:%d], D = [%d:%d:%d]\n",
            N, start_block_val, end_block_val, selection_block_step,
            start_double_block_val, end_double_block_val, selection_double_block_step);
    }

    double* left_matrix =  new double[N*N];
    double* right_matrix = new double[N*N];
    double* generator =    new double[N*N];
    double* left_copy =    new double[N*N];
    double* right_copy =   new double[N*N];

    const double memory_val = 40.0*N*N / (1024 * 1024);
    fprintf(log_file, "  > Объем выделенной памяти:   [ %.2lf Мб ]\n", memory_val);
    fprintf(log_file, "    ---------------------------------------------------------------\n");

    if (console_info_output)
    {
        printf("  > Объем выделенной памяти:   [ %.2lf Мб ]\n", memory_val);
        printf("    ---------------------------------------------------------------\n");
    }

    generate(left_matrix, N, N);
    generate(right_matrix, N, N);

    double min_res_time = DBL_MAX;
    double min_mul_time = DBL_MAX;
    double time_, realloc_time, mul_time, res_time;
    int optimal_block_by_wt = start_double_block_val;
    int optimal_block_by_mt = start_double_block_val;
    int optimal_d_block_by_wt = start_double_block_val;
    int optimal_d_block_by_mt = start_double_block_val;
        
    for (int block_size = start_block_val;
                      block_size <= end_block_val;
                      block_size += selection_block_step)
    {
        const int lower_bound = min(block_size - 1, start_double_block_val);
        const int upper_bound = min(block_size - 1, end_double_block_val);

        for (int double_block_size = lower_bound;
                          double_block_size <= upper_bound;
                          double_block_size += selection_double_block_step)
        {
            InitDispatchSystem();

            memcpy(left_copy,  left_matrix,  N*N*sizeof(double));
            memcpy(right_copy, right_matrix, N*N*sizeof(double));
            memset(generator,  0,            N*N*sizeof(double));

            TaskClass task_data = TaskClass(N, N,
                               block_size, block_size,
                               double_block_size, double_block_size);

            time_ = clock();
            standard_to_double_block_layout_reallocation(left_copy,
                N, N, block_size, block_size,
                double_block_size, double_block_size);
            standard_to_double_block_layout_reallocation(right_copy,
                N, N, block_size, block_size,
                double_block_size, double_block_size);
            realloc_time = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);

            time_ = clock();
            matrix_multiplication_double_block(generator,
                                               left_copy, right_copy,
                                               task_data, task_data);
            mul_time = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
            res_time = mul_time + realloc_time;

            fprintf(log_file, "  > B = %3d, D = %3d, MT = %f, RT = %f, WT = %f\n",
                    block_size, double_block_size, mul_time, realloc_time, res_time);

            if (console_info_output)
            {
                printf("  > B = %3d, D = %3d, MT = %f, RT = %f, WT = %f\n",
                    block_size, double_block_size, mul_time, realloc_time, res_time);
            }

            if (res_time < min_res_time)
            {
                min_res_time = res_time;
                optimal_block_by_wt = block_size;
                optimal_d_block_by_wt = double_block_size;
            }
            if (mul_time < min_mul_time)
            {
                min_mul_time = mul_time;
                optimal_block_by_mt = block_size;
                optimal_d_block_by_mt = double_block_size;
            }

            TurnOffDispatchSystem();
        }
    }
    fprintf(log_file, "  -----------------------------------------------------------------\n");
    fprintf(log_file, "  > Оптимальные параметры:\n");
    fprintf(log_file, "    По времени умножения (только умножение)          : %d и %d (%lf сек)\n",
            optimal_block_by_mt, optimal_d_block_by_mt, min_mul_time);
    fprintf(log_file, "    По общему времени (переразмещение + умножение)   : %d и %d (%lf сек)\n",
            optimal_block_by_wt, optimal_d_block_by_wt, min_res_time);

    fclose(log_file);

    if (console_info_output)
    {
        printf("  -----------------------------------------------------------------\n");
        printf("  > Оптимальные параметры:\n");
        printf("    По времени умножения (только умножение)          : %d и %d (%lf сек)\n",
            optimal_block_by_mt, optimal_d_block_by_mt, min_mul_time);
        printf("    По общему времени (переразмещение + умножение)   : %d и %d (%lf сек)\n",
            optimal_block_by_wt, optimal_d_block_by_wt, min_res_time);
    }

    delete[] left_copy;
    delete[] left_matrix;
    delete[] right_copy;
    delete[] right_matrix;
    delete[] generator;

    Blocks result;
    result.main_block = optimal_block_by_mt;
    result.small_block = optimal_d_block_by_mt;
    return result;
}

int select_optimal_block_size_multiplication(
                                       const int N,
                                       const int start_block_val,
                                       const int end_block_val,
                                       const int selection_step)
{
    InitDispatchSystem();

    double* left_matrix = new double[N*N];
    double* right_matrix = new double[N*N];
    double* generator = new double[N*N];
    double* left_copy = new double[N*N];
    double* right_copy = new double[N*N];

    const double memory_val = 40.0*N*N / (1024 * 1024);
    printf("\n [> Объем выделенной памяти:   [ %f Мб ]\n", memory_val);

    generate(left_matrix, N, N);
    generate(right_matrix, N, N);

    double min_res_time = DBL_MAX;
    double min_mul_time = DBL_MAX;
    double time_, realloc_time, mul_time, res_time;
    int optimal_by_whole_time = start_block_val;
    int optimal_by_mul_time = start_block_val;

    printf("\n [> Подбор размера блока...\n");
    printf("    Параметры: Размер матрицы = %d, перебор от %d до %d с шагом %d\n",
                  N, start_block_val, end_block_val, selection_step);
    for (int block_size = start_block_val; block_size <= end_block_val; block_size += selection_step)
    {
        memcpy(left_copy, left_matrix, N*N*sizeof(double));
        memcpy(right_copy, right_matrix, N*N*sizeof(double));
        
        TaskClass task_data = TaskClass(N, N, block_size, block_size);

        time_ = clock();
        standard_to_block_layout_reallocation(left_copy,
            N, N, block_size, block_size);
        standard_to_block_layout_reallocation(right_copy,
            N, N, block_size, block_size);
        realloc_time = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);

        time_ = clock();
        matrix_multiplication_block(generator, left_copy, right_copy,
                                    task_data, task_data);
        mul_time = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
        res_time = mul_time + realloc_time;        

        printf("  > B = %d, MT = %f сек, RT = %f сек, WT = %f сек\n",
                    block_size, mul_time, realloc_time, res_time);

        if (res_time < min_res_time)
        {
            min_res_time = res_time;
            optimal_by_whole_time = block_size;
        }
        if (mul_time < min_mul_time)
        {
            min_mul_time = mul_time;
            optimal_by_mul_time = block_size;
        }
    }
    printf("  -----------------------------------------------------------------\n");
    printf("  > OPTIMAL BLOCK SIZE: BY MT = %d, BY WT = %d\n\n", optimal_by_mul_time, optimal_by_whole_time);

    delete[] left_copy;
    delete[] left_matrix;
    delete[] right_copy;
    delete[] right_matrix;
    delete[] generator;

    TurnOffDispatchSystem();

    return optimal_by_mul_time;
}

void matrix_multiplication_tests(const TaskClass& params_left,
                          const TaskClass& params_right,
                          const bool       console_info_output)
{
    InitDispatchSystem();

    const int N1 = params_left.getDataRef().M_ROWS;
    const int N2 = params_left.getDataRef().M_COLS;
    const int N3 = params_right.getDataRef().M_COLS;
    const int B1 = params_left.getDataRef().B_ROWS;
    const int B2 = params_left.getDataRef().B_COLS;
    const int B3 = params_right.getDataRef().B_COLS;
    const int D1 = params_left.getDataRef().D_ROWS;
    const int D2 = params_left.getDataRef().D_COLS;
    const int D3 = params_right.getDataRef().D_COLS;

    FILE* log_file;
    fopen_s(&log_file, "../resource/log/mult_test_log.txt", "a");

    fprintf(log_file, "\n [> Запуск тестов для параметров:\n");
    fprintf(log_file, "    N1 = %5d, N2 = %5d, N3 = %5d\n", N1, N2, N3);
    fprintf(log_file, "    B1 = %5d, B2 = %5d, B3 = %5d\n", B1, B2, B3);
    fprintf(log_file, "    D1 = %5d, D2 = %5d, D3 = %5d\n", D1, D2, D3);

    if (console_info_output)
    {
        printf("\n\n [> Запуск тестов для параметров:\n");
        printf("    N1 = %5d, N2 = %5d, N3 = %5d\n", N1, N2, N3);
        printf("    B1 = %5d, B2 = %5d, B3 = %5d\n", B1, B2, B3);
        printf("    D1 = %5d, D2 = %5d, D3 = %5d\n", D1, D2, D3);
    }

    double* left_mat   = new double[N1*N2];
    double* rgt_mat    = new double[N2*N3];
    double* gen_mat_t  = new double[N1*N3];
    double* gen_mat_dt = new double[N1*N3];
    double* gen_mat_b  = new double[N1*N3];
    double* gen_mat_db = new double[N1*N3];

    const double memory_val = (8.0 * (N1*N2 + N2*N3 + 4*N1*N3)) / (1024 * 1024);
    fprintf(log_file, "  > Объем выделенной памяти:   [ %.2lf Мб ]\n", memory_val);
    fprintf(log_file, "    ---------------------------------------------------------------\n");
    if (console_info_output)
    {
        printf("  > Объем выделенной памяти:   [ %.2lf Мб ]\n", memory_val);
        printf("    ---------------------------------------------------------------\n");
    }

    fprintf(log_file, "  > Заполнение массивов...                         ");
    if (console_info_output)
    {
        printf("  > Заполнение массивов...                         ");
    }
    double time_ = clock();
    generate(left_mat, N1, N2);
    generate(rgt_mat, N2, N3);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Умножение матриц (тайлинг)...                  ");
    if (console_info_output)
    {
        printf("  > Умножение матриц (тайлинг)...                  ");
    }
    time_ = clock();
    matrix_multiplication_tiled(
        gen_mat_t, left_mat, rgt_mat,
        params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Умножение матриц (двойной тайлинг)...          ");
    if (console_info_output)
    {
        printf("  > Умножение матриц (двойной тайлинг)...          ");
    }
    time_ = clock();
    matrix_multiplication_double_tiled(
        gen_mat_dt, left_mat, rgt_mat,
        params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Блочное переразмещение...\n");
    fprintf(log_file, "    Левая матрица...                               ");
    if (console_info_output)
    {
        printf("  > Блочное переразмещение...\n");
        printf("    Левая матрица...                               ");
    }
    time_ = clock();
    standard_to_block_layout_reallocation(left_mat, N1, N2, B1, B2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "    Правая матрица...                              ");
    if (console_info_output)
    {
        printf("    Правая матрица...                              ");
    }
    time_ = clock();
    standard_to_block_layout_reallocation(rgt_mat, N2, N3, B2, B3);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Блочное умножение матриц...                    ");
    if (console_info_output)
    {
        printf("  > Блочное умножение матриц...                    ");
    }
    time_ = clock();
    matrix_multiplication_block(
        gen_mat_b, left_mat, rgt_mat,
        params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Обратное переразмещение...                     ");
    if (console_info_output)
    {
        printf("  > Обратное переразмещение...                     ");
    }
    time_ = clock();
    block_to_standard_layout_reallocation(gen_mat_b, N1, N3, B1, B3);
    block_to_standard_layout_reallocation(left_mat, N1, N2, B1, B2);
    block_to_standard_layout_reallocation(rgt_mat, N2, N3, B2, B3);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    TurnOffDispatchSystem();

    InitDispatchSystem();
    fprintf(log_file, "  > Двойное блочное переразмещение...\n");
    fprintf(log_file, "    Левая матрица...                               ");
    if (console_info_output)
    {
        printf("  > Двойное блочное переразмещение...\n");
        printf("    Левая матрица...                               ");
    }
    time_ = clock();
    standard_to_double_block_layout_reallocation(left_mat,
        N1, N2, B1, B2, D1, D2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "    Правая матрица...                              ");
    if (console_info_output)
    {
        printf("    Правая матрица...                              ");
    }
    time_ = clock();
    standard_to_double_block_layout_reallocation(rgt_mat,
        N2, N3, B2, B3, D2, D3);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Блочное умножение матриц (двойные блоки)...    ");
    if (console_info_output)
    {
        printf("  > Блочное умножение матриц (двойные блоки)...    ");
    }
    time_ = clock();
    matrix_multiplication_double_block(
        gen_mat_db, left_mat, rgt_mat,
        params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Обратное переразмещение...                     ");
    if (console_info_output)
    {
        printf("  > Обратное переразмещение...                     ");
    }
    time_ = clock();
    double_block_to_standard_layout_reallocation(gen_mat_db,
        N1, N3, B1, B3, D1, D3);
    double_block_to_standard_layout_reallocation(left_mat,
        N1, N2, B1, B2, D1, D2);
    double_block_to_standard_layout_reallocation(rgt_mat,
        N2, N3, B2, B3, D2, D3);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n\n", time_);
    }

    const double comparison_dbt  = compare_arrays(gen_mat_db, gen_mat_t,  N1*N3);
    const double comparison_dbdt = compare_arrays(gen_mat_db, gen_mat_dt, N1*N3);
    const double comparison_dbb  = compare_arrays(gen_mat_db, gen_mat_b,  N1*N3);
    const double comparison_bt   = compare_arrays(gen_mat_b,  gen_mat_t,  N1*N3);
    const double comparison_bdt  = compare_arrays(gen_mat_b,  gen_mat_dt, N1*N3);
    const double comparison_dtt  = compare_arrays(gen_mat_dt, gen_mat_t,  N1*N3);

    fprintf(log_file, "  > Попарное сравнение результатов:\n");
    fprintf(log_file, "    DIFFERENCE DB - T  [ %.6lf ]\n", comparison_dbt);
    fprintf(log_file, "    DIFFERENCE DB - DT [ %.6lf ]\n", comparison_dbdt);
    fprintf(log_file, "    DIFFERENCE DB - B  [ %.6lf ]\n", comparison_dbb);
    fprintf(log_file, "    DIFFERENCE B  - T  [ %.6lf ]\n", comparison_bt);
    fprintf(log_file, "    DIFFERENCE B  - DT [ %.6lf ]\n", comparison_bdt);
    fprintf(log_file, "    DIFFERENCE DT - T  [ %.6lf ]\n", comparison_dtt);
    if (console_info_output)
    {
        printf("  > Попарное сравнение результатов:\n");
        printf("    DIFFERENCE DB - T  [ %.6lf ]\n", comparison_dbt);
        printf("    DIFFERENCE DB - DT [ %.6lf ]\n", comparison_dbdt);
        printf("    DIFFERENCE DB - B  [ %.6lf ]\n", comparison_dbb);
        printf("    DIFFERENCE B  - T  [ %.6lf ]\n", comparison_bt);
        printf("    DIFFERENCE B  - DT [ %.6lf ]\n", comparison_bdt);
        printf("    DIFFERENCE DT - T  [ %.6lf ]\n", comparison_dtt);
    }
        
    fclose(log_file);

    delete[] left_mat;
    delete[] rgt_mat;
    delete[] gen_mat_t;
    delete[] gen_mat_dt;
    delete[] gen_mat_b;
    delete[] gen_mat_db;

    TurnOffDispatchSystem();
}

void floyd_test(const TaskClass& parameters,
                const bool console_info_output)
{
    InitDispatchSystem();

    const int N = parameters.getDataRef().M_ROWS;
    const int B = parameters.getDataRef().B_ROWS;

    FILE* log_file;
    fopen_s(&log_file, "../resource/log/floyd_test_log.txt", "a");

    fprintf(log_file, " [> Запуск тестов для параметров:\n");
    fprintf(log_file, "    N = %5d\n", N);
    fprintf(log_file, "    B = %5d\n", B);

    if (console_info_output)
    {
        printf(" [> Запуск тестов для параметров:\n");
        printf("    N = %5d\n", N);
        printf("    B = %5d\n", B);
    }

    double* mat_st = new double[N*N];
    double* mat_t  = new double[N*N];
    double* mat_b  = new double[N*N];
    
    const double memory_val = (24.0 * N*N) / (1024 * 1024);
    fprintf(log_file, "  > Объем выделенной памяти:   [ %.2lf Мб ]\n", memory_val);
    fprintf(log_file, "    ---------------------------------------------------------------\n");
    if (console_info_output)
    {
        printf("  > Объем выделенной памяти:   [ %.2lf Мб ]\n", memory_val);
        printf("    ---------------------------------------------------------------\n");
    }

    fprintf(log_file, "  > Заполнение массива...                          ");
    if (console_info_output)
    {
        printf("  > Заполнение массива...                          ");
    }
    double time_ = clock();
    generate(mat_st, N, N);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    memcpy(mat_t, mat_st, sizeof(double)*N*N);
    memcpy(mat_b, mat_st, sizeof(double)*N*N);

    fprintf(log_file, "  > Стандартный алгоритм Флойда...                 ");
    if (console_info_output)
    {
        printf("  > Стандартный алгоритм Флойда...                 ");
    }
    time_ = clock();
    floyd_standard(mat_st, N);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Алгоритм Флойда (тайлинг)...                   ");
    if (console_info_output)
    {
        printf("  > Алгоритм Флойда (тайлинг)...                   ");
    }
    time_ = clock();
    floyd_tiled(mat_t, N, B);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    
    fprintf(log_file, "  > Блочное переразмещение...                      ");
    if (console_info_output)
    {
        printf("  > Блочное переразмещение...                      ");
    }
    time_ = clock();
    standard_to_block_layout_reallocation(mat_b, N, N, B, B);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    
    fprintf(log_file, "  > Блочный алгоритм Флойда...                      ");
    if (console_info_output)
    {
        printf("  > Блочный алгоритм Флойда...                     ");
    }
    time_ = clock();
    floyd_block(mat_b, N, B);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Обратное переразмещение...                     ");
    if (console_info_output)
    {
        printf("  > Обратное переразмещение...                     ");
    }
    time_ = clock();
    block_to_standard_layout_reallocation(mat_b, N, N, B, B);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    
    const double comparison_bst = compare_arrays(mat_b, mat_st, N*N);
    const double comparison_bt  = compare_arrays(mat_b, mat_t, N*N);
    const double comparison_tst = compare_arrays(mat_t, mat_st, N*N);

    fprintf(log_file, "  > Попарное сравнение результатов:\n");
    fprintf(log_file, "    DIFFERENCE B - S  [ %.6lf ]\n", comparison_bst);
    fprintf(log_file, "    DIFFERENCE B - T  [ %.6lf ]\n", comparison_bt);
    fprintf(log_file, "    DIFFERENCE T - S  [ %.6lf ]\n\n", comparison_tst);
    if (console_info_output)
    {
        printf("  > Попарное сравнение результатов:\n");
        printf("    DIFFERENCE B - S  [ %.6lf ]\n", comparison_bst);
        printf("    DIFFERENCE B - T  [ %.6lf ]\n", comparison_bt);
        printf("    DIFFERENCE T - S  [ %.6lf ]\n\n", comparison_tst);
    }

    fclose(log_file);

    delete[] mat_st;
    delete[] mat_t;
    delete[] mat_b;

    TurnOffDispatchSystem();
}

int select_optimal_block_size_floyd(
    const int N,
    const int start_block_val,
    const int end_block_val,
    const int selection_step)
{
    double* original = new double[N*N];
    double* matrix = new double[N*N];
    
    const double memory_val = 24.0*N*N / (1024 * 1024);
    printf("\n [> Объем выделенной памяти:   [ %f Мб ]\n", memory_val);

    generate(original, N, N);

    double min_res_time = DBL_MAX;
    double min_alg_time = DBL_MAX;
    double time_, realloc_time, alg_time, res_time;
    int optimal_by_whole_time = start_block_val;
    int optimal_by_alg_time = start_block_val;

    printf("\n [> Подбор размера блока (алгоритм Флойда)...\n");
    printf("    Параметры: Размер матрицы = %d, перебор от %d до %d с шагом %d\n",
        N, start_block_val, end_block_val, selection_step);
    for (int block_size = start_block_val; block_size <= end_block_val; block_size += selection_step)
    {
        InitDispatchSystem();

        memcpy(matrix, original, N*N*sizeof(double));

        TaskClass task_data = TaskClass(N, N, block_size, block_size);

        time_ = clock();
        standard_to_block_layout_reallocation(matrix,
            N, N, block_size, block_size);
        realloc_time = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);

        time_ = clock();
        floyd_block(matrix, N, block_size);
        alg_time = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
        res_time = alg_time + realloc_time;

        printf("  > B = %d, AT = %f сек, RT = %f сек, WT = %f сек\n",
            block_size, alg_time, realloc_time, res_time);

        if (res_time < min_res_time)
        {
            min_res_time = res_time;
            optimal_by_whole_time = block_size;
        }
        if (alg_time < min_alg_time)
        {
            min_alg_time = alg_time;
            optimal_by_alg_time = block_size;
        }

        TurnOffDispatchSystem();
    }
    printf("  -----------------------------------------------------------------\n");
    printf("  > OPTIMAL BLOCK SIZE: BY MT = %d, BY WT = %d\n\n", optimal_by_alg_time, optimal_by_whole_time);

    delete[] matrix;
    delete[] original;

    return optimal_by_alg_time;
}

void reallocation_test(const TaskClass& parameters,
                       const bool console_info_output)
{
    InitDispatchSystem();

    const TaskData& data_ref = parameters.getDataRef();
    const int N1 = data_ref.M_ROWS;
    const int N2 = data_ref.M_COLS;
    const int B1 = data_ref.B_ROWS;
    const int B2 = data_ref.B_COLS;
    const int D1 = data_ref.D_ROWS;
    const int D2 = data_ref.D_COLS;

    FILE* log_file;
    fopen_s(&log_file, "../resource/log/reallocation_test_log.txt", "a");

    fprintf(log_file, "\n [> Запуск тестов для параметров:\n");
    fprintf(log_file, "    N1 = %5d, N2 = %5d\n", N1, N2);
    fprintf(log_file, "    B1 = %5d, B2 = %5d\n", B1, B2);
    fprintf(log_file, "    D1 = %5d, D2 = %5d\n", D1, D2);
    if (console_info_output)
    {
        printf("\n [> Запуск тестов для параметров:\n");
        printf("    N1 = %5d, N2 = %5d\n", N1, N2);
        printf("    B1 = %5d, B2 = %5d\n", B1, B2);
        printf("    D1 = %5d, D2 = %5d\n", D1, D2);
    }

    double* mat1 = new double[N1*N2];
    double* mat2 = new double[N1*N2];

    const double memory_val = (16.0 * N1*N2) / (1024 * 1024);
    fprintf(log_file, "  > Объем выделенной памяти:   [ %.2lf Мб ]\n", memory_val);
    fprintf(log_file, "    ---------------------------------------------------------------\n");
    if (console_info_output)
    {
        printf("  > Объем выделенной памяти:   [ %.2lf Мб ]\n", memory_val);
        printf("    ---------------------------------------------------------------\n");
    }

    fprintf(log_file, "  > Заполнение массива...                          ");
    if (console_info_output)
    {
        printf("  > Заполнение массива...                          ");
    }
    double time_ = clock();
    generate(mat1, N1, N2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    memcpy(mat2, mat1, N1*N2*sizeof(double));

    fprintf(log_file, "  > Блочное переразмещение...                      ");
    if (console_info_output)
    {
        printf("  > Блочное переразмещение...                      ");
    }
    time_ = clock();
    standard_to_block_layout_reallocation(mat1, N1, N2, B1, B2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    TurnOffDispatchSystem();

    InitDispatchSystem();
    fprintf(log_file, "  > Обратное переразмещение...                     ");
    if (console_info_output)
    {
        printf("  > Обратное переразмещение...                     ");
    }
    time_ = clock();
    block_to_standard_layout_reallocation(mat1, N1, N2, B1, B2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    TurnOffDispatchSystem();

    InitDispatchSystem();
    fprintf(log_file, "  > Двойное блочное переразмещение...              ");
    if (console_info_output)
    {
        printf("  > Двойное блочное переразмещение...              ");
    }
    time_ = clock();
    standard_to_double_block_layout_reallocation(mat1, N1, N2, B1, B2, D1, D2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    TurnOffDispatchSystem();

    InitDispatchSystem();
    fprintf(log_file, "  > Обратное переразмещение...                     ");
    if (console_info_output)
    {
        printf("  > Обратное переразмещение...                     ");
    }
    time_ = clock();
    double_block_to_standard_layout_reallocation(mat1, N1, N2, B1, B2, D1, D2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    TurnOffDispatchSystem();

    const double diff_norm = compare_arrays(mat1, mat2, N1*N2);
    fprintf(log_file, 
        "  > Норма разности:                                [ %.6lf ]\n\n",
        diff_norm);
    if (console_info_output)
    {
        printf("  > Норма разности:                                [ %.6lf ]\n\n",
            diff_norm);
    }

    fclose(log_file);

    delete[] mat1;
    delete[] mat2;
}

void qralg_test(const TaskClass& parameters,
    const bool console_info_output)
{
    InitDispatchSystem();

    const int N = parameters.getDataRef().M_ROWS;
    const int B1 = parameters.getDataRef().B_ROWS;
    const int B2 = parameters.getDataRef().B_COLS;
    const int D1 = parameters.getDataRef().D_ROWS;
    const int D2 = parameters.getDataRef().D_COLS;

    FILE* log_file;
    fopen_s(&log_file, "../resource/log/qralg_test_log.txt", "a");
    
    fprintf(log_file, "\n [> Запуск тестов для параметров:\n");
    fprintf(log_file, "    N = %5d\n", N);
    fprintf(log_file, "    B1 = %5d, B2 = %5d\n", B1, B2);
    fprintf(log_file, "    D1 = %5d, D2 = %5d\n", D1, D2);
    
    if (console_info_output)
    {
        printf("\n\n [> Запуск тестов для параметров:\n");
        printf("    N = %5d\n", N);
        printf("    B1 = %5d, B2 = %5d\n", B1, B2);
        printf("    D1 = %5d, D2 = %5d\n", D1, D2);
    }

    double* mat_t  = new double[N*N];
    double* mat_dt = new double[N*N];
    double* mat_b  = new double[N*N];
    double* mat_db = new double[N*N];

    const double memory_val = (32.0 * N * N) / (1024 * 1024);
    fprintf(log_file, "  > Объем выделенной памяти:   [ %.2lf Мб ]\n", memory_val);
    fprintf(log_file, "    ---------------------------------------------------------------\n");
    if (console_info_output)
    {
        printf("  > Объем выделенной памяти:   [ %.2lf Мб ]\n", memory_val);
        printf("    ---------------------------------------------------------------\n");
    }
    
    fprintf(log_file, "  > Заполнение массивов...                             ");
    if (console_info_output)
    {
        printf("  > Заполнение массивов...                             ");
    }
    double time_ = clock();
    generate(mat_t, N, N);
    memcpy(mat_dt, mat_t, N*N*sizeof(double));
    memcpy(mat_b,  mat_t, N*N*sizeof(double));
    memcpy(mat_db, mat_t, N*N*sizeof(double));
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > QR-алгоритм с тайлингом (WY-разложение)...         ");
    if (console_info_output)
    {
        printf("  > QR-алгоритм с тайлингом (WY-разложение)...         ");
    }
    time_ = clock();
    QR_WY_tiled(mat_t, parameters.getDataRef());
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > QR-алгоритм с двойным тайлингом (WY-разложение)... ");
    if (console_info_output)
    {
        printf("  > QR-алгоритм с двойным тайлингом (WY-разложение)... ");
    }
    time_ = clock();
    QR_WY_double_tiled(mat_dt, parameters.getDataRef());
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Блочное переразмещение...                          ");
    if (console_info_output)
    {
        printf("  > Блочное переразмещение...                          ");
    }
    time_ = clock();
    standard_to_block_layout_reallocation(mat_b, N, N, B1, B2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Блочный QR-алгоритм (WY-разложение)...             ");
    if (console_info_output)
    {
        printf("  > Блочный QR-алгоритм (WY-разложение)...             ");
    }
    time_ = clock();
    QR_WY_block(mat_b, parameters.getDataRef());
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Обратное переразмещение...                         ");
    if (console_info_output)
    {
        printf("  > Обратное переразмещение...                         ");
    }
    time_ = clock();
    block_to_standard_layout_reallocation(mat_b, N, N, B1, B2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    TurnOffDispatchSystem();
    
    InitDispatchSystem();
    fprintf(log_file, "  > Двойное блочное переразмещение...                  ");
    if (console_info_output)
    {
        printf("  > Двойное блочное переразмещение...                  ");
    }
    time_ = clock();
    standard_to_double_block_layout_reallocation(mat_db, N, N, B1, B2, D1, D2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    
    fprintf(log_file, "  > Двойной блочный QR-алгоритм (WY-разложение)...     ");
    if (console_info_output)
    {
        printf("  > Двойной блочный QR-алгоритм (WY-разложение)...     ");
    }
    time_ = clock();
    QR_WY_double_block(mat_db, parameters.getDataRef());
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n", time_);
    }
    
    fprintf(log_file, "  > Обратное переразмещение...                         ");
    if (console_info_output)
    {
        printf("  > Обратное переразмещение...                         ");
    }
    time_ = clock();
    double_block_to_standard_layout_reallocation(mat_db, N, N, B1, B2, D1, D2);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %.3lf секунд ]\n\n", time_);
    if (console_info_output)
    {
        printf("[ %.3lf секунд ]\n\n", time_);
    }

    const double comparison_dbt  = compare_arrays(mat_db, mat_t, N*N) ;
    const double comparison_dbdt = compare_arrays(mat_db, mat_dt, N*N);
    const double comparison_dbb  = compare_arrays(mat_db, mat_b, N*N);
    const double comparison_bt   = compare_arrays(mat_b, mat_t, N*N);
    const double comparison_bdt  = compare_arrays(mat_b, mat_dt, N*N);
    const double comparison_dtt  = compare_arrays(mat_dt, mat_t, N*N);

    fprintf(log_file, "  > Попарное сравнение результатов:\n");
    fprintf(log_file, "    DIFFERENCE DB - T  [ %.6lf ]\n", comparison_dbt);
    fprintf(log_file, "    DIFFERENCE DB - DT [ %.6lf ]\n", comparison_dbdt);
    fprintf(log_file, "    DIFFERENCE DB - B  [ %.6lf ]\n", comparison_dbb);
    fprintf(log_file, "    DIFFERENCE B  - T  [ %.6lf ]\n", comparison_bt);
    fprintf(log_file, "    DIFFERENCE B  - DT [ %.6lf ]\n", comparison_bdt);
    fprintf(log_file, "    DIFFERENCE DT - T  [ %.6lf ]\n", comparison_dtt);
    if (console_info_output)
    {
        printf("  > Попарное сравнение результатов:\n");
        printf("    DIFFERENCE DB - T  [ %.6lf ]\n", comparison_dbt);
        printf("    DIFFERENCE DB - DT [ %.6lf ]\n", comparison_dbdt);
        printf("    DIFFERENCE DB - B  [ %.6lf ]\n", comparison_dbb);
        printf("    DIFFERENCE B  - T  [ %.6lf ]\n", comparison_bt);
        printf("    DIFFERENCE B  - DT [ %.6lf ]\n", comparison_bdt);
        printf("    DIFFERENCE DT - T  [ %.6lf ]\n", comparison_dtt);
    }

    fclose(log_file);

    delete[] mat_t;
    delete[] mat_dt;
    delete[] mat_b;
    delete[] mat_db;

    TurnOffDispatchSystem();
}

