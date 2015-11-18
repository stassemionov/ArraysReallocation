#include <cstdlib>
#include <fstream>

#include "reallocating.h"
#include "service.h"

int main()
{
	setlocale(LC_ALL, "");
	ifstream file;
	file.open("../params.txt");
	// File needs to have a structure as follows (by a rows):
	// N1 - number of rows in matrix
	// N2 - number of columns in matrix
	// b1 - number of rows in block
	// b2 - number of columns in block
	int N1, N2, b1, b2;
	file >> N1 >> N2 >> b1 >> b2;
	file.close();

	cout << "N1 = " << N1 << endl;
	cout << "N2 = " << N2 << endl;
	cout << "b1 = " << b1 << endl;
	cout << "b2 = " << b2 << endl << endl;

	double* data = new double[N1*N2];	
	data = simple_fill(data, N1, N2);

	double* dt = get_reallocated(data, N1, N2, b1, b2);
	data = block_reallocate_matrix(data, N1, N2, b1, b2);

	cout << "Результат:\n";
//	print_to(cout, data, N1, N2, 4);
	cout << "Эталон:\n";
//	print_to(cout, dt,   N1, N2, 4);
	cout << "Норма разности: ";
	cout << compare_arrays(dt, data, N1*N2) << endl << endl;

	delete[] data;
	delete[] dt;
}