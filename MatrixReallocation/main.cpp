#include "reallocation.h"
#include "service.h"
#include "multiplication.h"
#include "testing.h"

#include <fstream>
#include <ctime>

#include "qralg.h"

using std::ifstream;

int main()
{
    setlocale(LC_ALL, "");
   
//    const TaskClass* floyd_params   = read_floyd_algorythm_parameters();
    const TaskClass* qr_params      = read_qr_parameters();
//    const TaskClass* mult_params    = read_multiplication_parameters();
//    const TaskClass* realloc_params = read_reallocation_test_parameters();

    const int N = qr_params->getDataRef().M_ROWS;
    const int B = qr_params->getDataRef().B_COLS;

    double* A = new double[N*N];
    double* Acpy = new double[N*N];
    double* Acpy2 = new double[N*N];

    generate(A, N, N);
    memcpy(Acpy, A, N*N*sizeof(double));
    memcpy(Acpy2, A, N*N*sizeof(double));

    printf("STANDARD START...\n"); 
    double time_ = clock();
    QR_WY_standard(Acpy, N, B);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    printf("STANDARD END: %lf секунд\n\n", time_);
 //   print_to(std::cout, Acpy, N, N, 10);

    printf("TILED START...\n");
    time_ = clock();
    QR_WY_tiled(A, qr_params->getDataRef());
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    printf("TILED END: %lf секунд\n\n", time_);
 //   print_to(std::cout, Acpy, N, N, 10);
    
    printf("CORRECTNESS: %lf\n\n", compare_arrays(A,Acpy,N*N));

    printf("REALLOCATION START...\n");
    time_ = clock();
    const BlockReallocationInfo* realloc_info = 
        standard_to_block_layout_reallocation(Acpy2, qr_params->getDataRef());
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    printf("REALLOCATION END: %lf секунд\n\n", time_);

    printf("BLOCK START...\n");
    time_ = clock();
    QR_WY_block(Acpy2, qr_params->getDataRef());
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    printf("BLOCK END: %lf секунд\n\n", time_);

    printf("INVERSE REALLOCATION START...\n");
    time_ = clock();
    block_to_standard_layout_reallocation(Acpy2, realloc_info);
    time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
    printf("INVERSE REALLOCATION END: %lf секунд\n\n", time_);
 //   print_to(std::cout, Acpy2, N, N, 10);

    printf("CORRECTNESS: %lf\n\n", compare_arrays(A, Acpy2, N*N));
}
