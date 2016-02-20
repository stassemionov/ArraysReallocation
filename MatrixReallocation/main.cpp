#include "reallocation.h"
#include "service.h"
#include "multiplication.h"
#include "testing.h"

#include <fstream>
#include <ctime>

using std::cout;
using std::ifstream;

int main()
{
    setlocale(LC_ALL, "");
    ostream& ostr = cout;
    ifstream file;
    file.open("../parameters.txt");
    // Task paramers file needs to have structure as follows:
    // <number of rows         in left       matrix>                  ['\n' | ' ']+
    // <number of columns/rows in left/rigth matrix>                  ['\n' | ' ']+
    // <number of columns      in right      matrix>                  ['\n' | ' ']+
    // <number of rows         in main block of left       matrix>    ['\n' | ' ']+
    // <number of columns/rows in main block of left/right matrix>    ['\n' | ' ']+
    // <number of columns      in main block of right      matrix>    ['\n' | ' ']+
    // <number of rows         in small block of left       matrix>   ['\n' | ' ']+
    // <number of columns/rows in small block of left/right matrix>   ['\n' | ' ']+
    // <number of rows         in small block of right      matrix>   ['\n' | ' ']+ <anything>
    int N1, N2, N3, b1, b2, b3, db1, db2, db3;
    file >> N1 >> N2 >> N3 >> b1 >> b2 >> b3 >> db1 >> db2 >> db3;
    file.close();

    ostr << "\n [> Параметры задачи:\n";
    ostr << "  > N1  = " << N1 << "\n";
    ostr << "  > N2  = " << N2 << "\n";
    ostr << "  > N3  = " << N3 << "\n";
    ostr << "  > b1  = " << b1 << "\n";
    ostr << "  > b2  = " << b2 << "\n";
    ostr << "  > b3  = " << b3 << "\n";
    ostr << "  > db1 = " << db1 << "\n";
    ostr << "  > db2 = " << db2 << "\n";
    ostr << "  > db3 = " << db3 << "\n\n";

    TaskClass params_left, params_right;
    params_left.makeData( N1, N2, b1, b2, db1, db2);
    params_right.makeData(N2, N3, b2, b3, db2, db3);

    matrix_multiplication_tests(params_left, params_right, true);

  /*  ostr << " [> Подбор параметров двойного блочного размещения...\n";
    Blocks optimal = select_optimal_double_block_size(1000,
                                                      15, 150, 1,
                                                      5, 40, 1,
                                                      false);
    ostr << "    Результат: " << optimal.main_block << " и " << optimal.small_block << std::endl;
    */

    
}
