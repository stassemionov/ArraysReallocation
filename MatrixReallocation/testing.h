/*

ostream& ostr = cout;
ifstream file;
file.open("../params.txt");
// Sizes data file needs to have structure as follows:
// <number of rows    in matrix>        ['\n' | ' ']+
// <number of columns in matrix>        ['\n' | ' ']+
// <number of rows    in main block>    ['\n' | ' ']+
// <number of columns in main block>    ['\n' | ' ']+
// <number of rows    in small block>   ['\n' | ' ']+
// <number of columns in small block>   [' ']+ <anything>
int N1, N2, b1, b2, db1, db2;
file >> N1 >> N2 >> b1 >> b2 >> db1 >> db2;
file.close();

ostr << "N1  = " << N1 << "\n";
ostr << "N2  = " << N2 << "\n";
ostr << "b1  = " << b1 << "\n";
ostr << "b2  = " << b2 << "\n";
ostr << "db1 = " << db1 << "\n";
ostr << "db2 = " << db2 << "\n\n";

TaskClass parameters;
parameters.makeData(N1, N2, b1, b2, db1, db2);

double* init_data = new double[N1*N2];
ostr << "Заполнение массива...     ";
double time_ = clock();
simple_fill(init_data, N1, N2);
time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
ostr << time_ << " секунд\n";
//  print_to(ostr, init_data, N1, N2, 4);

double* data = new double[N1*N2];
double* data_double = new double[N1*N2];
memcpy(data,        init_data, N1*N2*sizeof(double));
memcpy(data_double, init_data, N1*N2*sizeof(double));

ostr << "Вычисление эталона...     ";
time_ = clock();
double* standard_data = standard_to_block_layout_reallocation_buf(
data, parameters);
time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
ostr << time_ << " секунд\n";
//  print_to(ostr, standard_data, N1, N2, 4);

ostr << "Блочное переразмещение... ";
time_ = clock();
standard_to_block_layout_reallocation(data, parameters);
time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
ostr << time_ << " секунд\n\n";
//  print_to(ostr, data, N1, N2, 4);

ostr << "Норма разности: ";
ostr << compare_arrays(standard_data, data, N1*N2) << "\n\n";

ostr << "Вычисление эталона...     ";
time_ = clock();
double* standard_data_double =
standard_to_double_block_layout_reallocation_buf(data_double, parameters);
time_ = (clock() - time_) / (1.0 * CLOCKS_PER_SEC);
ostr << time_ << " секунд\n";
//  print_to(ostr, standard_data_double, N1, N2, 4);

ostr << "Двойное блочное переразмещение... ";
time_ = clock();
standard_to_double_block_layout_reallocation(data_double, parameters);
time_ = (clock() - time_) * 0.001;
ostr << time_ << " секунд\n\n";
//  print_to(ostr, data_double, N1, N2, 4);

ostr << "Норма разности: ";
ostr << compare_arrays(standard_data_double, data_double, N1*N2) << "\n\n";

*/