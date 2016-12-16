#ifndef _TASKDATA_H_
#define _TASKDATA_H_

#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

using std::vector;
using std::min;

// Stores main layout parameters
// like size of matrix and blocks,
// and some high-usage values,
// that are derived these parameters.
struct TaskData
{
    // * Main data
    int M_ROWS;          // matrix rows count
    int M_COLS;          // matrix columns count
    int B_ROWS;          // block rows count
    int B_COLS;          // block columns count
    int D_ROWS;          // double block rows count
    int D_COLS;          // double block columns count

    // * Additional data (derived from main data)
    int M_BLOCK_ROWS;    // number of block-stripes
    int M_BLOCK_COLS;    // number of block-columns
    int DIF_ROWS;        // number of rows in bottom-truncated block
    int DIF_COLS;        // number of columns in right-truncated block
    int MAIN_COLS;       // number of columns in matrix except truncated blocks
    int STRIPE_SIZE;     // number of elements in block-stripe
    int BLOCK_SIZE;      // number of elements in main block
    int DIF_BLOCK_SIZE;  // number of elements in right-truncated block
};

// Gives required interface over
// main layout parameters,
// that consists in index functions
// for transformations between
// standard and block layouts and inversely.
class TaskClass
{
public:

    TaskClass();

    TaskClass(const TaskData& init_data);

    TaskClass(const int& rows_count, const int& cols_count,
        const int& block_rows_count, const int& block_cols_count,
        const int& double_block_rows_count = 0,
        const int& double_block_cols_count = 0);

    const TaskData& getDataRef() const;

    // Принимает на вход адрес элемента в массиве со строчным размещением,
    // возвращает адрес этого элемента в массиве с двойным блочным размещением,
    // задаваемым набором параматров 'data'
    inline int _fastcall indexFunctionDbl(const int i_index,
                                const int j_index) const
    {
        // Координаты текущего большого блока
        const int&& b_row = i_index / m_data.B_ROWS;
        const int&& b_col = j_index / m_data.B_COLS;

        // Размеры текущего большого блока
        const int& bm = (b_row == m_data.M_BLOCK_ROWS - 1) ?
            ((m_data.DIF_ROWS != 0) ? m_data.DIF_ROWS : m_data.B_ROWS) : m_data.B_ROWS;
        const int& bn = (b_col == m_data.M_BLOCK_COLS - 1) ?
            ((m_data.DIF_COLS != 0) ? m_data.DIF_COLS : m_data.B_COLS) : m_data.B_COLS;

        // Размеры текущего малого блока
        const int& db1 = (bm < m_data.D_ROWS) ? bm : m_data.D_ROWS;
        const int& db2 = (bn < m_data.D_COLS) ? bn : m_data.D_COLS;

        // Количество малых блоков внутри текущего большого блока
        const int&& d_row_count = static_cast<int>(ceil(1.0 * bm / db1));
        const int&& d_col_count = static_cast<int>(ceil(1.0 * bn / db2));

        // Координаты элемента относительно текущего большого блока
        const int&& b_loc_i = i_index - b_row * m_data.B_ROWS;
        const int&& b_loc_j = j_index - b_col * m_data.B_COLS;

        // Координаты малого блока относительно большого, в котором он находится
        const int&& d_row = b_loc_i / db1;
        const int&& d_col = b_loc_j / db2;

        const int& current_small_block_h = (d_row == d_row_count - 1) ?
            ((bm % db1 != 0) ? bm % db1 : db1) : db1;
        const int& current_small_block_w = (d_col == d_col_count - 1) ?
            ((bn % db2 != 0) ? bn % db2 : db2) : db2;

        // Смещение начала текущего большого блока относительно начала массива:
        // b_shift = b_row * m_data.B_ROWS * m_data.M_COLS + bm * m_data.B_COLS * b_col

        // Смещение начала текущего малого блока
        // относительно начала текущего большого блока:
        // d_shift = d_row * db1 * bn + d_col * current_small_block_h * db2

        // Смещение элемента внутри текущего малого блока:
        // d_loc_i = b_loc_i - d_row * db1
        // d_loc_j = b_loc_j - d_col * db2
        // loc_shift = d_loc_i * current_small_block_w + d_loc_j

        // Адрес элемента в новом размещении складывается из
        // смещения большого блока (1), смещения малого блока (2) и
        // смещения внутри малого блока (3):
        // return b_shift + d_shift + loc_shift
        return
            b_row * m_data.B_ROWS * m_data.M_COLS + bm * m_data.B_COLS * b_col +  // 1

            d_row * db1 * bn + d_col * current_small_block_h * db2 +        // 2

            (b_loc_i - d_row * db1) * current_small_block_w +               // 3
            (b_loc_j - d_col * db2);
    }

    // Принимает на вход адрес элемента в массиве со строчным размещением,
    // возвращает адрес этого элемента в массиве с блочным размещением,
    // определяемым параметрами данного класса
    inline int _fastcall indexFunction(const int index) const
    {
        // Разложение адреса на строчную и столбцовую составляющие
        const int&& index_i = index / m_data.M_COLS;
        const int&& index_j = index % m_data.M_COLS;

        // Координаты текущего блока
        const int&& block_i = index_i / m_data.B_ROWS;
        const int&& block_j = index_j / m_data.B_COLS;

        // Размеры текущего блока
        const int& bm = (block_i == m_data.M_BLOCK_ROWS - 1) ?
            ((m_data.DIF_ROWS != 0) ? m_data.DIF_ROWS : m_data.B_ROWS) : m_data.B_ROWS;
        const int& bn = (block_j == m_data.M_BLOCK_COLS - 1) ?
            ((m_data.DIF_COLS != 0) ? m_data.DIF_COLS : m_data.B_COLS) : m_data.B_COLS;

        // Смещение текущего блока относительно начала массива:
        // int&& block_shift = block_i * m_data.B_ROWS * m_data.M_COLS +
        //                    block_j * bm * m_data.B_COLS;

        // Смещение элемента относительно начала текущего блока:
        // int&& loc_shift = bn * (index_i - block_i * m_data.B_ROWS) +
        //                  index_j - block_j * m_data.B_COLS;

        // Адрес элемента в новом размещении складывается из
        // смещения блока (1) и смещения элемента внутри блока (2):
        // return block_shift + loc_shift
        return block_i * m_data.B_ROWS * m_data.M_COLS +      // 1
            block_j * bm * m_data.B_COLS +                     // 1

            bn * (index_i - block_i * m_data.B_ROWS) +     // 2
            index_j - block_j * m_data.B_COLS;             // 2
    }

    // Принимает на вход адрес элемента в массиве со строчным размещением,
    // возвращает адрес этого элемента в массиве с блочным размещением,
    // определяемым параметрами данного класса,
    // при условии, что передан адрес элемента первой полосы блоков
    inline int _fastcall indexFunctionReduced(const int index) const
    {
        const int column = index % m_data.M_COLS;
        const int cur_block_width = (column < m_data.MAIN_COLS) ?
            m_data.B_COLS : m_data.DIF_COLS;
       
        // смещение в блочной полосе +
        // смещение в блоке этой полосы +
        // смещение в строке этого блока
        /*return m_data.BLOCK_SIZE * (column / m_data.B_COLS) +
            cur_block_width * (index / m_data.M_COLS) +
            column % m_data.B_COLS;*/
     
        // Optimized version.
        return m_data.B_ROWS * column +
            cur_block_width * (index / m_data.M_COLS) -
            (column % m_data.B_COLS) * (m_data.B_ROWS - 1);
    }

    inline int _fastcall indexFunctionReducedInverse(const int index) const
    {
        const int block_number = index / m_data.BLOCK_SIZE;
        const int local_shift = index % m_data.BLOCK_SIZE;
        const int cur_block_width =
            (m_data.STRIPE_SIZE - index - m_data.DIF_BLOCK_SIZE < 0) ?
            m_data.DIF_COLS : m_data.B_COLS;

        // смещение строками (от строки в блоке) +
        // смещение в строке (от числа блоков)
        // смещение в строке (от смещения в блоке)

        return (local_shift / cur_block_width) * m_data.M_COLS +
               block_number * m_data.B_COLS +
               local_shift % cur_block_width;
    }

private:
    TaskData m_data;
};

struct BlockReallocationInfo
{
    // Systems of Distinct Representatives of reallocation cycles
    vector<int> sdr_main;
    vector<int> sdr_addit;

    // Parameters of reallocation
    TaskClass main_data;
    TaskClass main_data_addit;

    ~BlockReallocationInfo() {}
};

struct DoubleBlockReallocationInfo
{
    const BlockReallocationInfo* upper_level_realloc_info = NULL;
    const BlockReallocationInfo* main_realloc_info = NULL;
    const BlockReallocationInfo* right_realloc_info = NULL;
    const BlockReallocationInfo* bottom_realloc_info = NULL;
    const BlockReallocationInfo* corner_realloc_info = NULL;
};

#endif  // _TASKDATA_H_
