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

    m_data.M_BLOCK_ROWS = static_cast<int>(ceil(1.0 * m_data.M_ROWS / m_data.B_ROWS));
    m_data.M_BLOCK_COLS = static_cast<int>(ceil(1.0 * m_data.M_COLS / m_data.B_COLS));
    m_data.DIF_COLS = m_data.M_COLS % m_data.B_COLS;
    m_data.DIF_ROWS = m_data.M_ROWS % m_data.B_ROWS;
    m_data.MAIN_COLS = m_data.M_COLS - m_data.DIF_COLS;
    m_data.STRIPE_SIZE = m_data.B_ROWS * m_data.M_COLS;
    m_data.BLOCK_SIZE = m_data.B_ROWS * m_data.B_COLS;
    m_data.DIF_BLOCK_SIZE = m_data.DIF_COLS * m_data.B_ROWS;
}

const TaskData& TaskClass::getDataRef() const
{
    return m_data;
}