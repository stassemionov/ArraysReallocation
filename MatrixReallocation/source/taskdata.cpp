#include "taskdata.h"

TaskClass::TaskClass() {}

TaskClass::TaskClass(const TaskData& init_data)
{
    TaskClass(init_data.M_ROWS, init_data.M_COLS,
        init_data.B_ROWS, init_data.B_COLS,
        init_data.D_ROWS, init_data.D_COLS);
}

TaskClass::TaskClass(const int& rows_count, const int& cols_count,
    const int& block_rows_count, const int& block_cols_count,
    const int& double_block_rows_count,
    const int& double_block_cols_count)
{
    m_data.M_ROWS = rows_count;
    m_data.M_COLS = cols_count;
    m_data.B_ROWS = block_rows_count;
    m_data.B_COLS = block_cols_count;
    m_data.D_ROWS = double_block_rows_count;
    m_data.D_COLS = double_block_cols_count;

    if (rows_count == 0       || cols_count == 0 ||
        block_rows_count == 0 || block_cols_count == 0)
    {
        return;
    }

    m_data.M_SIZE = m_data.M_ROWS * m_data.M_COLS;
    m_data.M_BLOCK_ROWS = static_cast<int>(ceil(1.0 * m_data.M_ROWS / m_data.B_ROWS));
    m_data.M_BLOCK_COLS = static_cast<int>(ceil(1.0 * m_data.M_COLS / m_data.B_COLS));
    m_data.DIF_COLS = m_data.M_COLS % m_data.B_COLS;
    m_data.DIF_ROWS = m_data.M_ROWS % m_data.B_ROWS;
    m_data.MAIN_COLS = m_data.M_COLS - m_data.DIF_COLS;
    m_data.STRIPE_SIZE = m_data.B_ROWS * m_data.M_COLS;
    m_data.COLUMN_SIZE = m_data.M_ROWS * m_data.B_COLS;
    m_data.BLOCK_SIZE = m_data.B_ROWS * m_data.B_COLS;
    m_data.DIF_ROWS_BLOCK_SIZE = m_data.DIF_ROWS * m_data.B_COLS;
    m_data.DIF_COLS_BLOCK_SIZE = m_data.DIF_COLS * m_data.B_ROWS;
    m_data.DIF_COLS_COLUMN_SIZE = m_data.DIF_COLS * m_data.M_ROWS;
}

const TaskData& TaskClass::getDataRef() const
{
    return m_data;
}

int TaskClass::indexFunctionDbl(const int i_index,
    const int j_index) const
{
    //  оординаты текущего большого блока
    const int b_row = i_index / m_data.B_ROWS;
    const int b_col = j_index / m_data.B_COLS;

    // –азмеры текущего большого блока
    const int bm = (b_row == m_data.M_BLOCK_ROWS - 1) ?
        ((m_data.DIF_ROWS != 0) ? m_data.DIF_ROWS : m_data.B_ROWS) : m_data.B_ROWS;
    const int bn = (b_col == m_data.M_BLOCK_COLS - 1) ?
        ((m_data.DIF_COLS != 0) ? m_data.DIF_COLS : m_data.B_COLS) : m_data.B_COLS;

    // –азмеры текущего малого блока
    const int db1 = (bm < m_data.D_ROWS) ? bm : m_data.D_ROWS;
    const int db2 = (bn < m_data.D_COLS) ? bn : m_data.D_COLS;

    //  оличество малых блоков внутри текущего большого блока
    const int d_row_count = static_cast<int>(ceil(1.0 * bm / db1));
    const int d_col_count = static_cast<int>(ceil(1.0 * bn / db2));

    //  оординаты элемента относительно текущего большого блока
    const int b_loc_i = i_index - b_row * m_data.B_ROWS;
    const int b_loc_j = j_index - b_col * m_data.B_COLS;

    //  оординаты малого блока относительно большого, в котором он находитс€
    const int d_row = b_loc_i / db1;
    const int d_col = b_loc_j / db2;

    const int current_small_block_h = (d_row == d_row_count - 1) ?
        ((bm % db1 != 0) ? bm % db1 : db1) : db1;
    const int current_small_block_w = (d_col == d_col_count - 1) ?
        ((bn % db2 != 0) ? bn % db2 : db2) : db2;

    // —мещение начала текущего большого блока относительно начала массива:
    // b_shift = b_row * m_data.B_ROWS * m_data.M_COLS + bm * m_data.B_COLS * b_col

    // —мещение начала текущего малого блока
    // относительно начала текущего большого блока:
    // d_shift = d_row * db1 * bn + d_col * current_small_block_h * db2

    // —мещение элемента внутри текущего малого блока:
    // d_loc_i = b_loc_i - d_row * db1
    // d_loc_j = b_loc_j - d_col * db2
    // loc_shift = d_loc_i * current_small_block_w + d_loc_j

    // јдрес элемента в новом размещении складываетс€ из
    // смещени€ большого блока (1), смещени€ малого блока (2) и
    // смещени€ внутри малого блока (3):
    // return b_shift + d_shift + loc_shift
    return
        b_row * m_data.B_ROWS * m_data.M_COLS + bm * m_data.B_COLS * b_col +  // 1

        d_row * db1 * bn + d_col * current_small_block_h * db2 +        // 2

        (b_loc_i - d_row * db1) * current_small_block_w +               // 3
        (b_loc_j - d_col * db2);
}

int TaskClass::indexFunction(const int index) const
{
    // –азложение адреса на строчную и столбцовую составл€ющие
    const int index_i = index / m_data.M_COLS;
    const int index_j = index % m_data.M_COLS;

    //  оординаты текущего блока
    const int block_i = index_i / m_data.B_ROWS;
    const int block_j = index_j / m_data.B_COLS;

    // –азмеры текущего блока
    const int bm = (block_i == m_data.M_BLOCK_ROWS - 1) ?
        ((m_data.DIF_ROWS != 0) ? m_data.DIF_ROWS : m_data.B_ROWS) : m_data.B_ROWS;
    const int bn = (block_j == m_data.M_BLOCK_COLS - 1) ?
        ((m_data.DIF_COLS != 0) ? m_data.DIF_COLS : m_data.B_COLS) : m_data.B_COLS;

    // —мещение текущего блока относительно начала массива:
    // int&& block_shift = block_i * m_data.B_ROWS * m_data.M_COLS +
    //                    block_j * bm * m_data.B_COLS;

    // —мещение элемента относительно начала текущего блока:
    // int&& loc_shift = bn * (index_i - block_i * m_data.B_ROWS) +
    //                  index_j - block_j * m_data.B_COLS;

    // јдрес элемента в новом размещении складываетс€ из
    // смещени€ блока (1) и смещени€ элемента внутри блока (2):
    // return block_shift + loc_shift
    return block_i * m_data.B_ROWS * m_data.M_COLS +      // 1
        block_j * bm * m_data.B_COLS +                     // 1

        bn * (index_i - block_i * m_data.B_ROWS) +     // 2
        index_j - block_j * m_data.B_COLS;             // 2
}

int TaskClass::indexFunctionTransposed(const int index) const
{
    const int row = index / m_data.M_COLS;
    const int col = index % m_data.M_COLS;
    const int brow = row / m_data.B_ROWS;
    const int bcol = col / m_data.B_COLS;
    const int cur_block_width = min(m_data.B_COLS, m_data.M_COLS - bcol * m_data.B_COLS);
    const int cur_block_height = min(m_data.B_ROWS, m_data.M_ROWS - brow * m_data.B_ROWS);
    return bcol * m_data.COLUMN_SIZE +
        (row / m_data.B_ROWS) * (m_data.B_ROWS * cur_block_width) +
        (col % m_data.B_COLS) * cur_block_height +
        row % m_data.B_ROWS;
}

int TaskClass::indexFunctionTransposedInverse(const int index) const
{
    // Shift of block-columns.
    const int stripe_shift = index / m_data.COLUMN_SIZE;
    // Width of current block.
    const int cur_block_width =
        min(m_data.B_COLS, m_data.M_COLS - stripe_shift*m_data.B_COLS);
    const int loc_index_block_column =
        index - stripe_shift * m_data.COLUMN_SIZE;
    // Shift of blocks in current block-column.
    const int block_col_shift =
        loc_index_block_column / (m_data.B_ROWS * cur_block_width);
    const int cur_block_height =
        min(m_data.B_ROWS, m_data.M_ROWS - block_col_shift*m_data.B_ROWS);
    const int loc_index_block =
        loc_index_block_column - cur_block_width*m_data.B_ROWS * block_col_shift;
    const int loc_hor_shift = loc_index_block / cur_block_height;
    const int loc_index_col = loc_index_block - loc_hor_shift * cur_block_height;

    return (block_col_shift * m_data.B_ROWS + loc_index_col) * m_data.M_COLS +
        stripe_shift * m_data.B_COLS + loc_hor_shift;
}