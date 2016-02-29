#include "reallocation.h"
#include "service.h"
#include "multiplication.h"
#include "testing.h"

#include <fstream>
#include <ctime>

using std::ifstream;

int main()
{
    setlocale(LC_ALL, "");
   
    const TaskClass* floyd_params   = read_floyd_algorythm_parameters();
    const TaskClass* mult_params    = read_multiplication_parameters();
    const TaskClass* realloc_params = read_reallocation_test_parameters();


}
