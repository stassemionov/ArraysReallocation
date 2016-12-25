#ifndef _TASKDATA_H_
#define _TASKDATA_H_

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
    int DIF_ROWS_BLOCK_SIZE;  // number of elements in bottom-truncated block
    int DIF_COLS_BLOCK_SIZE;  // number of elements in right-truncated block
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
    int _fastcall indexFunctionDbl(const int i_index,
        const int j_index) const;

    // Принимает на вход адрес элемента в массиве со строчным размещением,
    // возвращает адрес этого элемента в массиве с блочным размещением,
    // определяемым параметрами данного класса
    int _fastcall indexFunction(const int index) const;

    // Принимает на вход адрес элемента в массиве со строчным размещением,
    // возвращает адрес этого элемента в массиве с блочным размещением,
    // определяемым параметрами данного класса,
    // при условии, что передан адрес элемента первой полосы блоков
    inline int indexFunctionReduced(const int index) const
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

    inline int indexFunctionReducedInverse(const int index) const
    {
        const int block_number = index / m_data.BLOCK_SIZE;
        const int local_shift = index % m_data.BLOCK_SIZE;
        const int cur_block_width =
            (m_data.STRIPE_SIZE - index - m_data.DIF_COLS_BLOCK_SIZE < 0) ?
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
    
    inline bool isTrivial() const
    {
        return sdr_main.empty() || (main_data.getDataRef().BLOCK_SIZE == 0);
    }
};

struct DoubleBlockReallocationInfo
{
    const BlockReallocationInfo* upper_level_realloc_info = nullptr;
    const BlockReallocationInfo* main_realloc_info        = nullptr;
    const BlockReallocationInfo* right_realloc_info       = nullptr;
    const BlockReallocationInfo* bottom_realloc_info      = nullptr;
    const BlockReallocationInfo* corner_realloc_info      = nullptr;

    // This part is unused when transposition option is disable.
    const BlockReallocationInfo* main_transpose_info   = nullptr;
    const BlockReallocationInfo* right_transpose_info  = nullptr;
    const BlockReallocationInfo* bottom_transpose_info = nullptr;
    const BlockReallocationInfo* corner_transpose_info = nullptr;
};

#endif  // _TASKDATA_H_
