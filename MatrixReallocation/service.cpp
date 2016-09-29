#include "service.h"

#include <iomanip>
#include <fstream>
#include <algorithm>

using std::min;
using std::max;
using std::setw;
using std::ifstream;
using std::make_pair;

double* generate(double* data_ptr, 
    const int rows_count, const int cols_count,
    const double lbound, const double ubound)
{
    double rnd = (ubound - lbound) / RAND_MAX;
    srand(static_cast<int>(time(NULL)));
    const size_t length = static_cast<size_t>(rows_count * cols_count);
    for (size_t i = 0; i < length; ++i)
    {
        data_ptr[i] = lbound + rand() * rnd;
    }
    return data_ptr;
}

double* simple_fill(double* data_ptr,
    const int rows_count, const int cols_count)
{
    size_t i(0);
    const size_t length = static_cast<size_t>(rows_count * cols_count);
    while (i < length)
    {
        data_ptr[i] = double(++i);
    }
    return data_ptr;
}

void print_to(ostream& ostr, const double* data_ptr,
    const int rows_count, const int cols_count, const int place)
{
    const double* row = data_ptr;
    for (int i = 0; i < rows_count; ++i)
    {
        for (int j = 0; j < cols_count; ++j)
        {
            ostr << setw(place) << row[j];
        }
        row += cols_count;
        ostr << '\n';
    }
    ostr << '\n';
}

double _fastcall compare_arrays(
    const double* data1, const double* data2, const size_t length)
{
    double s(0.0);
    for (size_t i = 0; i < length; ++i)
    {
        s += fabs(data1[i] - data2[i]);
    }
    return s;
}

int _fastcall gcd(const int u, const int v)
{
    // simple cases (termination)
    if (u == v)
        return u;

    if (u == 0)
        return v;

    if (v == 0)
        return u;

    // look for factors of 2
    if (~u & 1)     // u is even
    {
        if (v & 1)  // v is odd
            return gcd(u >> 1, v);
        else    // both u and v are even
            return gcd(u >> 1, v >> 1) << 1;
    }

    if (~v & 1)     // u is odd, v is even
        return gcd(u, v >> 1);

    // reduce larger argument
    if (u > v)
        return gcd((u - v) >> 1, v);

    return gcd((v - u) >> 1, u);
}

bool _fastcall bin_search(const int val, const vector<int>& vec)
{
    if (vec.empty())
    {
        return false;
    }
    if (val > vec.back())
    {
        return false;
    }
    if (val < vec.front())
    {
        return false;
    }

    int l(0);
    int r = static_cast<int>(vec.size()) - 1;
    int m(0);
    while (l < r)
    {
        m = (l + r) >> 1;
        if (vec[m] < val)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    return vec[r] == val;
}

pair<TaskClass, TaskClass> read_multiplication_parameters(
    const string& file_name)
{
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

    ifstream file;
    file.open(file_name);

    int N1, N2, N3, b1, b2, b3, db1, db2, db3;
    file >> N1 >> N2 >> N3 >> b1 >> b2 >> b3 >> db1 >> db2 >> db3;
    file.close();

    return make_pair<TaskClass, TaskClass>(
        TaskClass(N1, N2, b1, b2, db1, db2),
        TaskClass(N2, N3, b2, b3, db2, db3));
}

TaskClass read_reallocation_test_parameters(const string& file_name)
{
    // Task paramers file needs to have structure as follows:
    // <number of rows    in matrix>                  ['\n' | ' ']+
    // <number of columns in matrix>                  ['\n' | ' ']+
    // <number of rows    in main block of matrix>    ['\n' | ' ']+
    // <number of columns in main block of matrix>    ['\n' | ' ']+
    // <number of rows    in small block of matrix>   ['\n' | ' ']+
    // <number of rows    in small block of matrix>   ['\n' | ' ']+ <anything>

    ifstream file;
    file.open(file_name);

    int N1, N2, b1, b2, db1, db2;
    file >> N1 >> N2 >> b1 >> b2 >> db1 >> db2;
    file.close();

    return TaskClass(N1, N2, b1, b2, db1, db2);
}

TaskClass read_floyd_algorythm_parameters(const string& file_name)
{
    // Task paramers file needs to have structure as follows:
    // <number of rows/columns in matrix>   ['\n' | ' ']+
    // <number of rows in block of matrix>  ['\n' | ' ']+ <anything>
    
    ifstream file;
    file.open(file_name);

    int N, b;
    file >> N >> b;
    file.close();

    return TaskClass(N, N, b, b);
}

TaskClass read_qr_parameters(const string& file_name)
{
    // Task paramers file needs to have structure as follows:
    // <number of rows/columns in matrix>                  ['\n' | ' ']+
    // <number of rows         in main block of matrix>    ['\n' | ' ']+
    // <number of columns      in main block of matrix>    ['\n' | ' ']+
    // <number of rows         in small block of matrix>   ['\n' | ' ']+
    // <number of columns      in small block of matrix>   ['\n' | ' ']+ <anything>

    ifstream file;
    file.open(file_name);

    int N, b1, b2, db1, db2;
    file >> N >> b1 >> b2 >> db1 >> db2;
    file.close();

    return TaskClass(N, N, b1, b2, db1, db2);
}

void copy_minor_double_block(double* dst_ptr,
    const double* src_ptr,
    const TaskData& layout_data,
    const int row_coord,
    const int col_coord)
{
    // Main parameters of data layout
    const int N = layout_data.M_ROWS;
    const int b1 = layout_data.B_ROWS;
    const int b2 = layout_data.B_COLS;
    const int db1 = layout_data.D_ROWS;
    const int db2 = layout_data.D_COLS;

    if (row_coord >= N || col_coord >= N)
    {
        return;
    }

    // Big blocks counts in matrix
    const int block_count_in_col = static_cast<int>(ceil(1.0 * N / b1));
    const int block_count_in_row = static_cast<int>(ceil(1.0 * N / b2));
    // Small blocks counts in big blocks
    const int dblock_count_in_bcol = static_cast<int>(ceil(1.0 * b1 / db1));
    const int dblock_count_in_brow = static_cast<int>(ceil(1.0 * b2 / db2));
    const int dblock_count_in_diff_bcol = (N % b1 == 0) ?
        dblock_count_in_bcol :
        static_cast<int>(ceil(1.0 * (N % b1) / db1));
    const int dblock_count_in_diff_brow = (N % b2 == 0) ?
        dblock_count_in_brow :
        static_cast<int>(ceil(1.0 * (N % b2) / db2));

    // Indices of big block, that contains element (row_coord, col_coord)
    const int curr_i_block_ind = row_coord / b1;
    const int curr_j_block_ind = col_coord / b2;

    // Horizontal and vertical shifts to corner big block
    const int d1_shift = curr_i_block_ind * b1;
    const int d2_shift = curr_j_block_ind * b2;
    // Real values of height/width of current block-stripe/column
    const int real_b1 = min(b1, N - d1_shift);
    const int real_b2 = min(b2, N - d2_shift);
    // Counts of small blocks in corner big block
    const int dblock_count_in_corner_block_in_col =
        (curr_i_block_ind == block_count_in_col - 1) ?
        dblock_count_in_diff_bcol : dblock_count_in_bcol;
    const int dblock_count_in_corner_block_in_row =
        (curr_j_block_ind == block_count_in_row - 1) ?
        dblock_count_in_diff_brow : dblock_count_in_brow;

    // Corner coordinates in its big block
    const int loc_row_coord = row_coord - d1_shift;
    const int loc_col_coord = col_coord - d2_shift;
    // Local coordinates of small block containing element (row_coord, col_coord)
    const int id_border = loc_row_coord / db1;
    const int jd_border = loc_col_coord / db2;
    // Local shifts to small block containing element (row_coord, col_coord)
    // and its sizes
    const int id0_shift = id_border * db1;
    const int id0_block_h = min(db1, real_b1 - id0_shift);
    const int jd0_shift = jd_border * db2;
    const int jd0_block_w = min(db2, real_b2 - jd0_shift);
    // Corner coordinates in its small block
    const int dloc_row_coord = loc_row_coord - id0_shift;
    const int dloc_col_coord = loc_col_coord - jd0_shift;

    double* dst_stripe = dst_ptr + d1_shift * N;
    const double* src_stripe = src_ptr + (dst_stripe - dst_ptr);

    // Copying of corner incomplete big block
    {
        double* dst_block = dst_stripe + d2_shift * real_b1;
        const double* src_block = src_stripe + (dst_block - dst_stripe);
        // Pointers to small-block-stripe containing element (row_coord, col_coord)
        double* dst0_dstripe = dst_block + id0_shift * real_b2;
        const double* src0_dstripe = src_block + (dst0_dstripe - dst_block);
        // Firstly, copy incomplete small-block-stripe
        for (int jd = jd_border; jd < dblock_count_in_corner_block_in_row; ++jd)
        {
            const int j_bd_shift = jd * db2;
            const int j_lb_loc = max(loc_col_coord, j_bd_shift) - j_bd_shift;
            const int jd_block_w = min(db2, real_b2 - j_bd_shift);

            double* dst_dblock = dst0_dstripe + j_bd_shift * id0_block_h;
            const double* src_dblock = src0_dstripe + (dst_dblock - dst0_dstripe);
            
            const int cpy_size = (jd_block_w - j_lb_loc) * sizeof(double);
            for (int l = dloc_row_coord; l < id0_block_h; ++l)
            {
                const int shift = l*jd_block_w + j_lb_loc;
                memcpy(dst_dblock + shift,
                    src_dblock + shift,
                    cpy_size);
            }
        }

        for (int id = id_border + 1; id < dblock_count_in_corner_block_in_col; ++id)
        {
            const int id_shift = id * db1;
            const int id_block_h = min(db1, real_b1 - id_shift);
            double* dst_dstripe = dst_block + id_shift * real_b2;
            const double* src_dstripe = src_block + (dst_dstripe - dst_block);

            double* dst_dblock = dst_dstripe + jd0_shift * id_block_h;
            const double* src_dblock = src_dstripe + (dst_dblock - dst_dstripe);

            const int cpy_size = (jd0_block_w - dloc_col_coord) * sizeof(double);
            for (int l = 0; l < id_block_h; ++l)
            {
                const int shift = l*jd0_block_w + dloc_col_coord;
                memcpy(dst_dblock + shift,
                    src_dblock + shift,
                    cpy_size);
            }

            const int shift = (jd0_shift + jd0_block_w) * id_block_h;
            memcpy(dst_dstripe + shift,
                src_dstripe + shift,
                (real_b2*id_block_h - shift)*sizeof(double));
        }
    }

    // Copy passing over first big-block-stripe, which can be incomplete
    for (int jb = curr_j_block_ind + 1; jb < block_count_in_row; ++jb)
    {
        const int j_mb_shift = jb * b2;
        const int jb_block_w = min(b2, N - j_mb_shift);
        double* dst_block = dst_stripe + j_mb_shift * real_b1;
        const double* src_block = src_stripe + (dst_block - dst_stripe);

        double* dst0_dstripe = dst_block + id0_shift * jb_block_w;
        const double* src0_dstripe = src_block + (dst0_dstripe - dst_block);

        const int jd_ub = (jb == block_count_in_row - 1) ?
            dblock_count_in_diff_brow : dblock_count_in_brow;
        for (int jd = 0; jd < jd_ub; ++jd)
        {
            const int j_bd_shift = jd * db2;
            const int jd_block_w = min(db2, jb_block_w - j_bd_shift);
            double* dst_dblock = dst0_dstripe + j_bd_shift * id0_block_h;
            const double* src_dblock = src0_dstripe + (dst_dblock - dst0_dstripe);

            const int shift = dloc_row_coord * jd_block_w;
            const int cpy_size = (id0_block_h*jd_block_w - shift)*sizeof(double);
            memcpy(dst_dblock + shift,
                   src_dblock + shift,
                   cpy_size);
        }

        const int shift = (id0_shift + id0_block_h) * jb_block_w;
        const int cpy_size = (real_b1*jb_block_w - shift)*sizeof(double);
        memcpy(dst_block + shift, src_block + shift, cpy_size);
    }

    // Copy passing over left big-block-column, which can be incomplete
    for (int ib = curr_i_block_ind + 1; ib < block_count_in_col; ++ib)
    {
        const int i_mb_shift = ib * b1;
        const int ib_block_h = min(b1, N - i_mb_shift);
        double* dst_block = dst_ptr + i_mb_shift*N + d2_shift*ib_block_h;
        const double* src_block = src_ptr + (dst_block - dst_ptr);
        
        const int id_ub = (ib == block_count_in_col - 1) ?
            dblock_count_in_diff_bcol : dblock_count_in_bcol;
        for (int id = 0; id < id_ub; ++id)
        {
            const int i_bd_shift = id * db1;
            const int id_block_h = min(db1, ib_block_h - i_bd_shift);
            double* dst_dstripe = dst_block + i_bd_shift * real_b2;
            const double* src_dstripe = src_block + (dst_dstripe - dst_block);

            double* dst_dblock_shifted = dst_dstripe + jd0_shift*id_block_h + dloc_col_coord;
            const double* src_dblock_shifted = src_dstripe + (dst_dblock_shifted - dst_dstripe);
            const int cpy_line_size = (jd0_block_w - dloc_col_coord) * sizeof(double);
            for (int l = 0; l < id_block_h; ++l)
            {
                const int shift = l * jd0_block_w;
                memcpy(dst_dblock_shifted + shift,
                       src_dblock_shifted + shift,
                       cpy_line_size);
            }

            const int shift = id_block_h * (jd0_shift + jd0_block_w);
            const int cpy_size = (id_block_h*real_b2 - shift)*sizeof(double);
            memcpy(dst_dstripe + shift,
                   src_dstripe + shift,
                   cpy_size);
        }
    }

    // Copy passing over every big-block-stripe, which is complete
    for (int ib = curr_i_block_ind + 1; ib < block_count_in_col; ++ib)
    {
        const int i_mb_shift = ib * b1;
        const int ib_block_h = min(b1, N - i_mb_shift);
        const int shift = ib_block_h*(d2_shift + real_b2);

        double* dst_block = dst_ptr + i_mb_shift*N + shift;
        const double* src_block = src_ptr + (dst_block - dst_ptr);
        const int cpy_size = (ib_block_h*N - shift)*sizeof(double);
        memcpy(dst_block, src_block, cpy_size);
    }
}

void set_minor_double_block(double* dst_ptr,
    const int val,
    const TaskData& layout_data,
    const int row_coord,
    const int col_coord)
{
    // Main parameters of data layout
    const int N = layout_data.M_ROWS;
    const int b1 = layout_data.B_ROWS;
    const int b2 = layout_data.B_COLS;
    const int db1 = layout_data.D_ROWS;
    const int db2 = layout_data.D_COLS;

    if (row_coord >= N || col_coord >= N)
    {
        return;
    }

    // Big blocks counts in matrix
    const int block_count_in_col = static_cast<int>(ceil(1.0 * N / b1));
    const int block_count_in_row = static_cast<int>(ceil(1.0 * N / b2));
    // Small blocks counts in big blocks
    const int dblock_count_in_bcol = static_cast<int>(ceil(1.0 * b1 / db1));
    const int dblock_count_in_brow = static_cast<int>(ceil(1.0 * b2 / db2));
    const int dblock_count_in_diff_bcol = (N % b1 == 0) ?
        dblock_count_in_bcol :
        static_cast<int>(ceil(1.0 * (N % b1) / db1));
    const int dblock_count_in_diff_brow = (N % b2 == 0) ?
        dblock_count_in_brow :
        static_cast<int>(ceil(1.0 * (N % b2) / db2));

    // Indices of big block, that contains element (row_coord, col_coord)
    const int curr_i_block_ind = row_coord / b1;
    const int curr_j_block_ind = col_coord / b2;

    // Horizontal and vertical shifts to corner big block
    const int d1_shift = curr_i_block_ind * b1;
    const int d2_shift = curr_j_block_ind * b2;
    // Real values of height/width of current block-stripe/column
    const int real_b1 = min(b1, N - d1_shift);
    const int real_b2 = min(b2, N - d2_shift);
    // Counts of small blocks in corner big block
    const int dblock_count_in_corner_block_in_col =
        (curr_i_block_ind == block_count_in_col - 1) ?
        dblock_count_in_diff_bcol : dblock_count_in_bcol;
    const int dblock_count_in_corner_block_in_row =
        (curr_j_block_ind == block_count_in_row - 1) ?
        dblock_count_in_diff_brow : dblock_count_in_brow;

    // Corner coordinates in its big block
    const int loc_row_coord = row_coord - d1_shift;
    const int loc_col_coord = col_coord - d2_shift;
    // Local coordinates of small block containing element (row_coord, col_coord)
    const int id_border = loc_row_coord / db1;
    const int jd_border = loc_col_coord / db2;
    // Local shifts to small block containing element (row_coord, col_coord)
    // and its sizes
    const int id0_shift = id_border * db1;
    const int id0_block_h = min(db1, real_b1 - id0_shift);
    const int jd0_shift = jd_border * db2;
    const int jd0_block_w = min(db2, real_b2 - jd0_shift);
    // Corner coordinates in its small block
    const int dloc_row_coord = loc_row_coord - id0_shift;
    const int dloc_col_coord = loc_col_coord - jd0_shift;

    double* dst_stripe = dst_ptr + d1_shift * N;

    // Copying of corner incomplete big block
    {
        double* dst_block = dst_stripe + d2_shift * real_b1;
        // Pointers to small-block-stripe containing element (row_coord, col_coord)
        double* dst0_dstripe = dst_block + id0_shift * real_b2;
        // Firstly, copy incomplete small-block-stripe
        for (int jd = jd_border; jd < dblock_count_in_corner_block_in_row; ++jd)
        {
            const int j_bd_shift = jd * db2;
            const int j_lb_loc = max(loc_col_coord, j_bd_shift) - j_bd_shift;
            const int jd_block_w = min(db2, real_b2 - j_bd_shift);
            double* dst_dblock_shifted = dst0_dstripe + j_bd_shift * id0_block_h + j_lb_loc;

            const int cpy_size = (jd_block_w - j_lb_loc) * sizeof(double);
            for (int l = dloc_row_coord; l < id0_block_h; ++l)
            {
                memset(dst_dblock_shifted + l*jd_block_w, val, cpy_size);
            }
        }

        for (int id = id_border + 1; id < dblock_count_in_corner_block_in_col; ++id)
        {
            const int id_shift = id * db1;
            const int id_block_h = min(db1, real_b1 - id_shift);
            double* dst_dstripe = dst_block + id_shift * real_b2;
            double* dst_dblock_shifted = dst_dstripe + jd0_shift * id_block_h + dloc_col_coord;

            const int cpy_size = (jd0_block_w - dloc_col_coord) * sizeof(double);
            for (int l = 0; l < id_block_h; ++l)
            {
                memset(dst_dblock_shifted + l*jd0_block_w, val, cpy_size);
            }

            const int shift = (jd0_shift + jd0_block_w) * id_block_h;
            memset(dst_dstripe + shift, val,
                (real_b2*id_block_h - shift)*sizeof(double));
        }
    }

    // Copy passing over first big-block-stripe, which can be incomplete
    for (int jb = curr_j_block_ind + 1; jb < block_count_in_row; ++jb)
    {
        const int j_mb_shift = jb * b2;
        const int jb_block_w = min(b2, N - j_mb_shift);
        double* dst_block = dst_stripe + j_mb_shift * real_b1;
        double* dst0_dstripe = dst_block + id0_shift * jb_block_w;

        const int jd_ub = (jb == block_count_in_row - 1) ?
            dblock_count_in_diff_brow : dblock_count_in_brow;
        for (int jd = 0; jd < jd_ub; ++jd)
        {
            const int j_bd_shift = jd * db2;
            const int jd_block_w = min(db2, jb_block_w - j_bd_shift);
            double* dst_dblock = dst0_dstripe + j_bd_shift * id0_block_h;

            const int shift = dloc_row_coord * jd_block_w;
            memset(dst_dblock + shift, val,
                (id0_block_h*jd_block_w - shift)*sizeof(double));
        }

        const int shift = (id0_shift + id0_block_h) * jb_block_w;
        const int set_size = (real_b1*jb_block_w - shift)*sizeof(double);
        memset(dst_block + shift, val, set_size);
    }

    // Copy passing over left big-block-column, which can be incomplete
    for (int ib = curr_i_block_ind + 1; ib < block_count_in_col; ++ib)
    {
        const int i_mb_shift = ib * b1;
        const int ib_block_h = min(b1, N - i_mb_shift);
        double* dst_block = dst_ptr + i_mb_shift*N + d2_shift*ib_block_h;

        const int id_ub = (ib == block_count_in_col - 1) ?
            dblock_count_in_diff_bcol : dblock_count_in_bcol;
        for (int id = 0; id < id_ub; ++id)
        {
            const int i_bd_shift = id * db1;
            const int id_block_h = min(db1, ib_block_h - i_bd_shift);
            double* dst_dstripe = dst_block + i_bd_shift * real_b2;

            double* dst_dblock_shifted = dst_dstripe + jd0_shift*id_block_h + dloc_col_coord;
            const int set_line_size = (jd0_block_w - dloc_col_coord) * sizeof(double);
            for (int l = 0; l < id_block_h; ++l)
            {
                memset(dst_dblock_shifted + l*jd0_block_w, val, set_line_size);
            }

            const int shift = id_block_h * (jd0_shift + jd0_block_w);
            memset(dst_dstripe + shift, val,
                (id_block_h*real_b2 - shift)*sizeof(double));
        }
    }

    // Copy passing over every big-block-stripe, which is complete
    for (int ib = curr_i_block_ind + 1; ib < block_count_in_col; ++ib)
    {
        const int i_mb_shift = ib * b1;
        const int ib_block_h = min(b1, N - i_mb_shift);
        const int shift = ib_block_h*(d2_shift + real_b2);
        memset(dst_ptr + i_mb_shift*N + shift, val, (ib_block_h*N - shift)*sizeof(double));
    }
}

void copy_minor_block(double* dst_ptr,
    const double* src_ptr,
    const TaskData& layout_data,
    const int row_coord,
    const int col_coord)
{
    // Main parameters of data layout
    const int N = layout_data.M_ROWS;
    const int b1 = layout_data.B_ROWS;
    const int b2 = layout_data.B_COLS;

    if (row_coord >= N || col_coord >= N)
    {
        return;
    }

    // Big blocks counts in matrix
    const int block_count_in_col = static_cast<int>(ceil(1.0 * N / b1));
    const int block_count_in_row = static_cast<int>(ceil(1.0 * N / b2));
    // Indices of big block, that contains element (row_coord, col_coord)
    const int curr_i_block_ind = row_coord / b1;
    const int curr_j_block_ind = col_coord / b2;
    // Horizontal and vertical shifts to corner big block
    const int d1_shift = curr_i_block_ind * b1;
    const int d2_shift = curr_j_block_ind * b2;
    // Real values of height/width of current block-stripe/column
    const int real_b1 = min(b1, N - d1_shift);
    const int real_b2 = min(b2, N - d2_shift);
    // Corner coordinates in its big block
    const int loc_row_coord = row_coord - d1_shift;
    const int loc_col_coord = col_coord - d2_shift;

    double* dst_stripe = dst_ptr + d1_shift * N;
    const double* src_stripe = src_ptr + (dst_stripe - dst_ptr);

    // * Copy of corner incomplete block

    // Pointers to first element to copy in corner block
    const double* src_row_shifted_CORNER = src_stripe + real_b1*d2_shift +
        loc_row_coord*real_b2 + loc_col_coord;
    double* dst_row_shifted_CORNER = dst_stripe +
        (src_row_shifted_CORNER - src_stripe);
    const size_t copy_size_CORNER = (real_b2 - loc_col_coord) * sizeof(double);
    for (int i = loc_row_coord; i < real_b1; ++i)
    {
        memcpy(dst_row_shifted_CORNER,
            src_row_shifted_CORNER,
            copy_size_CORNER);
        dst_row_shifted_CORNER += real_b2;
        src_row_shifted_CORNER += real_b2;
    }

    // * Copy of first big-block-stripe, which can be incomplete
    for (int jb = curr_j_block_ind + 1; jb < block_count_in_row; ++jb)
    {
        const int jb_shift = jb * b2;
        const int jb_block_width = min(b2, N - jb_shift);
        const size_t copy_size = jb_block_width * sizeof(double);
        const double* src_row = src_stripe + jb_shift * real_b1 +
            loc_row_coord * jb_block_width;
        double* dst_row = dst_stripe + (src_row - src_stripe);
        for (int i = loc_row_coord; i < real_b1; ++i)
        {
            memcpy(dst_row, src_row, copy_size);
            dst_row += jb_block_width;
            src_row += jb_block_width;
        }
    }

    // * Copy of first big-block-column, which can be incomplete
    const size_t copy_size_LEFT = (real_b2 - loc_col_coord) * sizeof(double);
    for (int ib = curr_i_block_ind + 1; ib < block_count_in_col; ++ib)
    {
        const int ib_shift = ib * b1;
        const int ib_block_h = min(b1, N - ib_shift);
        double* dst_row = dst_ptr + ib_shift * N +
            ib_block_h * d2_shift + loc_col_coord;
        const double* src_row = src_ptr + (dst_row - dst_ptr);

        for (int i = 0; i < ib_block_h; ++i)
        {
            memcpy(dst_row, src_row, copy_size_LEFT);
            dst_row += real_b2;
            src_row += real_b2;
        }
    }

    // * Copy of remaining part of source matrix
    const int main_part_shift = d2_shift + real_b2;
    const int main_part_width_in_bytes = (N - main_part_shift) * sizeof(double);
    for (int ib = curr_i_block_ind + 1; ib < block_count_in_col; ++ib)
    {
        const int ib_shift = ib * b1;
        const int ib_block_h = min(b1, N - ib_shift);
        double* dst_block = dst_ptr + ib_shift * N +
            ib_block_h * main_part_shift;
        const double* src_block = src_ptr + (dst_block - dst_ptr);

        memcpy(dst_block, src_block, ib_block_h * main_part_width_in_bytes);
    }
}

void set_minor_block(double* dst_ptr,
    const int val,
    const TaskData& layout_data,
    const int row_coord,
    const int col_coord)
{
    // Main parameters of data layout
    const int N = layout_data.M_ROWS;
    const int b1 = layout_data.B_ROWS;
    const int b2 = layout_data.B_COLS;

    if (row_coord >= N || col_coord >= N)
    {
        return;
    }

    // Big blocks counts in matrix
    const int block_count_in_col = static_cast<int>(ceil(1.0 * N / b1));
    const int block_count_in_row = static_cast<int>(ceil(1.0 * N / b2));
    // Indices of big block, that contains element (row_coord, col_coord)
    const int curr_i_block_ind = row_coord / b1;
    const int curr_j_block_ind = col_coord / b2;
    // Horizontal and vertical shifts to corner big block
    const int d1_shift = curr_i_block_ind * b1;
    const int d2_shift = curr_j_block_ind * b2;
    // Real values of height/width of current block-stripe/column
    const int real_b1 = min(b1, N - d1_shift);
    const int real_b2 = min(b2, N - d2_shift);
    // Corner coordinates in its big block
    const int loc_row_coord = row_coord - d1_shift;
    const int loc_col_coord = col_coord - d2_shift;

    double* dst_stripe = dst_ptr + d1_shift * N;

    // * Set value for elements of corner incomplete block

    // Pointer to first element to copy in corner block
    double* dst_row_shifted_CORNER = dst_stripe + real_b1*d2_shift +
        loc_row_coord*real_b2 + loc_col_coord;
    const size_t copy_size_CORNER = (real_b2 - loc_col_coord) * sizeof(double);
    for (int i = loc_row_coord; i < real_b1; ++i)
    {
        memset(dst_row_shifted_CORNER, val, copy_size_CORNER);
        dst_row_shifted_CORNER += real_b2;
    }

    // * Set value for elements of first big-block-stripe, which can be incomplete
    for (int jb = curr_j_block_ind + 1; jb < block_count_in_row; ++jb)
    {
        const int jb_shift = jb * b2;
        const int jb_block_width = min(b2, N - jb_shift);
        const size_t copy_size = jb_block_width * sizeof(double);
        double* dst_row = dst_stripe + jb_shift * real_b1 +
            loc_row_coord * jb_block_width;
        for (int i = loc_row_coord; i < real_b1; ++i)
        {
            memset(dst_row, 0, copy_size);
            dst_row += jb_block_width;
        }
    }

    // * Copy of first big-block-column, which can be incomplete
    const size_t copy_size_LEFT = (real_b2 - loc_col_coord) * sizeof(double);
    for (int ib = curr_i_block_ind + 1; ib < block_count_in_col; ++ib)
    {
        const int ib_shift = ib * b1;
        const int ib_block_h = min(b1, N - ib_shift);
        double* dst_row = dst_ptr + ib_shift * N +
            ib_block_h * d2_shift + loc_col_coord;

        for (int i = 0; i < ib_block_h; ++i)
        {
            memset(dst_row, 0, copy_size_LEFT);
            dst_row += real_b2;
        }
    }

    // * Copy of remaining part of source matrix
    const int main_part_shift = d2_shift + real_b2;
    const int main_part_width_in_bytes = (N - main_part_shift) * sizeof(double);
    for (int ib = curr_i_block_ind + 1; ib < block_count_in_col; ++ib)
    {
        const int ib_shift = ib * b1;
        const int ib_block_h = min(b1, N - ib_shift);
        double* dst_block = dst_ptr + ib_shift * N +
            ib_block_h * main_part_shift;

        memset(dst_block, 0, ib_block_h * main_part_width_in_bytes);
    }
}