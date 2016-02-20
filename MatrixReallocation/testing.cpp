#include "testing.h"
#include "reallocation.h"
#include "multiplication.h"
#include "service.h"

#include <cstring>
#include <ctime>
#include <algorithm>

using std::min;

Blocks select_optimal_double_block_size(const unsigned int N,
    const unsigned int start_block_val,
    const unsigned int end_block_val,
    const unsigned int selection_block_step,
    const unsigned int start_double_block_val,
    const unsigned int end_double_block_val,
    const unsigned int selection_double_block_step,
    const bool         console_info_output)
{
    FILE* log_file;
    fopen_s(&log_file, "../double_block_selection_test_log.txt", "w+");

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
    unsigned int optimal_block_by_wt = start_double_block_val;
    unsigned int optimal_block_by_mt = start_double_block_val;
    unsigned int optimal_d_block_by_wt = start_double_block_val;
    unsigned int optimal_d_block_by_mt = start_double_block_val;
        
    for (unsigned int block_size = start_block_val;
                      block_size <= end_block_val;
                      block_size += selection_block_step)
    {
        const unsigned int lower_bound = min(block_size - 1, start_double_block_val);
        const unsigned int upper_bound = min(block_size - 1, end_double_block_val);

        for (unsigned int double_block_size = lower_bound;
                          double_block_size <= upper_bound;
                          double_block_size += selection_double_block_step)
        {
            memcpy(left_copy,  left_matrix,  N*N*sizeof(double));
            memcpy(right_copy, right_matrix, N*N*sizeof(double));
            memset(generator,  0,            N*N*sizeof(double));

            TaskClass task_data;
            task_data.makeData(N, N,
                               block_size, block_size,
                               double_block_size, double_block_size);

            time_ = clock();
            standard_to_double_block_layout_reallocation(left_copy, task_data);
            standard_to_double_block_layout_reallocation(right_copy, task_data);
            realloc_time = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);

            time_ = clock();
            block_matrix_multiplication_double_tiled(generator, 
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

unsigned int select_optimal_block_size(const unsigned int N,
                                       const unsigned int start_block_val,
                                       const unsigned int end_block_val,
                                       const unsigned int selection_step)
{
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
    unsigned int optimal_by_whole_time = start_block_val;
    unsigned int optimal_by_mul_time = start_block_val;

    printf("\n [> Подбор размера блока...\n");
    printf("    Параметры: Размер матрицы = %d, перебор от %d до %d с шагом %d\n",
                  N, start_block_val, end_block_val, selection_step);
    for (unsigned int block_size = start_block_val; block_size <= end_block_val; block_size += selection_step)
    {
        memcpy(left_copy, left_matrix, N*N*sizeof(double));
        memcpy(right_copy, right_matrix, N*N*sizeof(double));
        
        TaskClass task_data;
        task_data.makeData(N, N, block_size, block_size);

        time_ = clock();
        standard_to_block_layout_reallocation(left_copy, task_data);
        standard_to_block_layout_reallocation(right_copy, task_data);
        realloc_time = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);

        time_ = clock();
        block_matrix_multiplication_tiled(generator, left_copy, right_copy,
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

    return optimal_by_mul_time;
}

void matrix_multiplication_tests(const TaskClass& params_left,
                          const TaskClass& params_right,
                          const bool       console_info_output)
{
    const int N1 = params_left.getDataRef().M_ROWS;
    const int N2 = params_left.getDataRef().M_COLS;
    const int N3 = params_right.getDataRef().M_COLS;

    FILE* log_file;
    fopen_s(&log_file, "../mult_test_log.txt", "w");

    fprintf(log_file, "\n\n [> Запуск тестов для параметров:\n");
    fprintf(log_file, "    N1 = %5d, N2 = %5d, N3 = %5d\n", N1, N2, N3);
    fprintf(log_file, "    B1 = %3d, B2 = %3d, B3 = %3d\n",
        params_left.getDataRef().B_ROWS, params_left.getDataRef().B_COLS, params_right.getDataRef().B_COLS);
    fprintf(log_file, "    D1 = %3d, D2 = %3d, D3 = %3d\n",
        params_left.getDataRef().D_ROWS, params_left.getDataRef().D_COLS, params_right.getDataRef().D_COLS);

    if (console_info_output)
    {
        printf("\n\n [> Запуск тестов для параметров:\n");
        printf("    N1 = %5d, N2 = %5d, N3 = %5d\n", N1, N2, N3);
        printf("    B1 = %3d, B2 = %3d, B3 = %3d\n",
           params_left.getDataRef().B_ROWS, params_left.getDataRef().B_COLS, params_right.getDataRef().B_COLS);
        printf("    D1 = %3d, D2 = %3d, D3 = %3d\n",
            params_left.getDataRef().D_ROWS, params_left.getDataRef().D_COLS, params_right.getDataRef().D_COLS);
    }

    double* left_mat = new double[N1*N2];
    double* rgt_mat = new double[N2*N3];
    double* left_mat_copy = new double[N1*N2];
    double* rgt_mat_copy = new double[N2*N3];
    double* gen_mat = new double[N1*N3];

    const double memory_val = (8.0 * (2*N1*N2 + 2*N2*N3 + N1*N3)) / (1024 * 1024);
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
    fprintf(log_file, "[ %lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %lf секунд ]\n", time_);
    }

    fprintf(log_file, "  > Умножение матриц (тайлинг)...                  ");
    if (console_info_output)
    {
        printf("  > Умножение матриц (тайлинг)...                  ");
    }
    time_ = clock();
    matrix_multiplication_tiled(gen_mat, left_mat, rgt_mat,
        params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %lf секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %lf секунд ]\n", time_);
    }

    memset(gen_mat, 0, N1*N3*sizeof(double));

    fprintf(log_file, "  > Умножение матриц (двойной тайлинг)...          ");
    if (console_info_output)
    {
        printf("  > Умножение матриц (двойной тайлинг)...          ");
    }
    time_ = clock();
    matrix_multiplication_double_tiled(gen_mat, left_mat, rgt_mat,
                                       params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %lf  секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %lf  секунд ]\n", time_);
    }

    memset(gen_mat, 0, N1*N3*sizeof(double));
    
    memcpy(left_mat_copy, left_mat, N1*N2*sizeof(double));
    fprintf(log_file, "  > Блочное переразмещение...\n");
    fprintf(log_file, "    Левая матрица...                               ");
    if (console_info_output)
    {
        printf("  > Блочное переразмещение...\n");
        printf("    Левая матрица...                               ");
    }
    time_ = clock();
    standard_to_block_layout_reallocation(left_mat, params_left);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %lf  секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %lf  секунд ]\n", time_);
    }

    memcpy(rgt_mat_copy, rgt_mat, N2*N3*sizeof(double));
    fprintf(log_file, "    Правая матрица...                              ");
    if (console_info_output)
    {
        printf("    Правая матрица...                              ");
    }
    time_ = clock();
    standard_to_block_layout_reallocation(rgt_mat, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %lf  секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %lf  секунд ]\n", time_);
    }

    fprintf(log_file, "  > Блочное умножение матриц (тайлинг)...          ");
    if (console_info_output)
    {
        printf("  > Блочное умножение матриц (тайлинг)...          ");
    }
    time_ = clock();
    block_matrix_multiplication_tiled(gen_mat, left_mat, rgt_mat,
                                      params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %lf  секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %lf  секунд ]\n", time_);
    }

    fprintf(log_file, "  > Двойное блочное переразмещение...\n");
    fprintf(log_file, "    Левая матрица...                               ");
    if (console_info_output)
    {
        printf("  > Двойное блочное переразмещение...\n");
        printf("    Левая матрица...                               ");
    }
    time_ = clock();
    standard_to_double_block_layout_reallocation(left_mat_copy, params_left);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %lf  секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %lf  секунд ]\n", time_);
    }

    fprintf(log_file, "    Правая матрица...                              ");
    if (console_info_output)
    {
        printf("    Правая матрица...                              ");
    }
    time_ = clock();
    standard_to_double_block_layout_reallocation(rgt_mat_copy, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %lf  секунд ]\n", time_);
    if (console_info_output)
    {
        printf("[ %lf  секунд ]\n", time_);
    }

    fprintf(log_file, "  > Блочное умножение матриц (двойной тайлинг)...  ");
    if (console_info_output)
    {
        printf("  > Блочное умножение матриц (двойной тайлинг)...  ");
    }
    time_ = clock();
    block_matrix_multiplication_double_tiled(gen_mat, left_mat_copy, rgt_mat_copy,
                                             params_left, params_right);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    fprintf(log_file, "[ %lf  секунд ]\n\n", time_);
    if (console_info_output)
    {
        printf("[ %lf  секунд ]\n\n", time_);
    }
    
    fclose(log_file);
}

void correctness_test(const TaskClass& params_left,
    const TaskClass& params_right,
    const bool       console_info_output = false)
{

}