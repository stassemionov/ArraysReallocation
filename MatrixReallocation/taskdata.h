#ifndef _TASKDATA_H_
#define _TASKDATA_H_

#include <cstdlib>
#include <cmath>

struct TaskData
{
    // Main data
    int M_ROWS;     // matrix rows count
    int M_COLS;     // matrix columns count
    int B_ROWS;     // block rows count
    int B_COLS;     // block columns count
    int D_ROWS;     // double block rows count
    int D_COLS;     // double block columns count

    // Additional data (derived from main data)
    int M_BLOCK_ROWS;       // ����� ������ � ����������� �������
    int M_BLOCK_COLS;       // ����� ������ � ����������� ������
    int DIF_ROWS;           // �������� ����� �� �������
    int DIF_COLS;           // �������� ����� �� ������
};

class TaskClass
{
public:

    TaskClass();

    TaskClass(const TaskData& init_data);

    void makeData(const int& rows_count, const int& cols_count,
        const int& block_rows_count, const int& block_cols_count,
        const int& double_block_rows_count = 0,
        const int& double_block_cols_count = 0);

    const TaskData& getDataRef() const;

    // ��������� �� ���� ����� �������� � ������� �� �������� �����������,
    // ���������� ����� ����� �������� � ������� � ������� ������� �����������,
    // ���������� ������� ���������� 'data'
    inline int indexFunctionDbl(const int& i_index,
                                const int& j_index) const
    {
        // ���������� �������� �������� �����
        const int&& b_row = i_index / m_data.B_ROWS;
        const int&& b_col = j_index / m_data.B_COLS;

        // ������� �������� �������� �����
        const int& bm = (b_row == m_data.M_BLOCK_ROWS - 1) ?
            ((m_data.DIF_ROWS != 0) ? m_data.DIF_ROWS : m_data.B_ROWS) : m_data.B_ROWS;
        const int& bn = (b_col == m_data.M_BLOCK_COLS - 1) ?
            ((m_data.DIF_COLS != 0) ? m_data.DIF_COLS : m_data.B_COLS) : m_data.B_COLS;

        // ������� �������� ������ �����
        const int& db1 = (bm < m_data.D_ROWS) ? bm : m_data.D_ROWS;
        const int& db2 = (bn < m_data.D_COLS) ? bn : m_data.D_COLS;

        // ���������� ����� ������ ������ �������� �������� �����
        const int&& d_row_count = static_cast<int>(ceil(1.0 * bm / db1));
        const int&& d_col_count = static_cast<int>(ceil(1.0 * bn / db2));

        // ���������� �������� ������������ �������� �������� �����
        const int&& b_loc_i = i_index - b_row * m_data.B_ROWS;
        const int&& b_loc_j = j_index - b_col * m_data.B_COLS;

        // ���������� ������ ����� ������������ ��������, � ������� �� ���������
        const int&& d_row = b_loc_i / db1;
        const int&& d_col = b_loc_j / db2;

        const int& current_small_block_h = (d_row == d_row_count - 1) ?
            ((bm % db1 != 0) ? bm % db1 : db1) : db1;
        const int& current_small_block_w = (d_col == d_col_count - 1) ?
            ((bn % db2 != 0) ? bn % db2 : db2) : db2;

        // �������� ������ �������� �������� ����� ������������ ������ �������:
        // b_shift = b_row * m_data.B_ROWS * m_data.M_COLS + bm * m_data.B_COLS * b_col

        // �������� ������ �������� ������ �����
        // ������������ ������ �������� �������� �����:
        // d_shift = d_row * db1 * bn + d_col * current_small_block_h * db2

        // �������� �������� ������ �������� ������ �����:
        // d_loc_i = b_loc_i - d_row * db1
        // d_loc_j = b_loc_j - d_col * db2
        // loc_shift = d_loc_i * current_small_block_w + d_loc_j

        // ����� �������� � ����� ���������� ������������ ��
        // �������� �������� ����� (1), �������� ������ ����� (2) �
        // �������� ������ ������ ����� (3):
        // return b_shift + d_shift + loc_shift
        return
            b_row * m_data.B_ROWS * m_data.M_COLS + bm * m_data.B_COLS * b_col +  // 1

            d_row * db1 * bn + d_col * current_small_block_h * db2 +        // 2

            (b_loc_i - d_row * db1) * current_small_block_w +               // 3
            (b_loc_j - d_col * db2);
    }

    // ��������� �� ���� ����� �������� � ������� �� �������� �����������,
    // ���������� ����� ����� �������� � ������� � ������� �����������,
    // ���������� ������� ���������� 'data'
    inline int indexFunction(const int& index) const
    {
        // ���������� ������ �� �������� � ���������� ������������
        const int&& index_i = index / m_data.M_COLS;
        const int&& index_j = index % m_data.M_COLS;

        // ���������� �������� �����
        const int&& block_i = index_i / m_data.B_ROWS;
        const int&& block_j = index_j / m_data.B_COLS;

        // ������� �������� �����
        const int& bm = (block_i == m_data.M_BLOCK_ROWS - 1) ?
            ((m_data.DIF_ROWS != 0) ? m_data.DIF_ROWS : m_data.B_ROWS) : m_data.B_ROWS;
        const int& bn = (block_j == m_data.M_BLOCK_COLS - 1) ?
            ((m_data.DIF_COLS != 0) ? m_data.DIF_COLS : m_data.B_COLS) : m_data.B_COLS;

        // �������� �������� ����� ������������ ������ �������:
        // int&& block_shift = block_i * m_data.B_ROWS * m_data.M_COLS +
        //                    block_j * bm * m_data.B_COLS;

        // �������� �������� ������������ ������ �������� �����:
        // int&& loc_shift = bn * (index_i - block_i * m_data.B_ROWS) +
        //                  index_j - block_j * m_data.B_COLS;

        // ����� �������� � ����� ���������� ������������ ��
        // �������� ����� (1) � �������� �������� ������ ����� (2):
        // return block_shift + loc_shift
        return block_i * m_data.B_ROWS * m_data.M_COLS +      // 1
            block_j * bm * m_data.B_COLS +                     // 1

            bn * (index_i - block_i * m_data.B_ROWS) +     // 2
            index_j - block_j * m_data.B_COLS;             // 2
    }

    // ��������� �� ���� ����� �������� � ������� �� �������� �����������,
    // ���������� ����� ����� �������� � ������� � ������� �����������,
    // ���������� ������� ���������� 'data',
    // ��� �������, ��� ������� ����� �������� ������ ������ ������
    inline int indexFunctionReduced(const int& index) const
    {
        const int&& i_col = index % m_data.M_COLS;
      //  const int& curr_block_width = 
      //      (i_col < m_data.M_COLS - m_data.DIF_COLS) ?
      //      m_data.B_COLS : m_data.DIF_COLS;

        if (i_col < m_data.M_COLS - m_data.DIF_COLS)
        {
            return m_data.B_ROWS * m_data.B_COLS * (i_col / m_data.B_COLS) +
                m_data.B_COLS * (index / m_data.M_COLS) +
                i_col % m_data.B_COLS;
        }
        else
        {
            return m_data.B_ROWS * m_data.B_COLS * (i_col / m_data.B_COLS) +
                m_data.DIF_COLS * (index / m_data.M_COLS) +
                i_col % m_data.B_COLS;
        }
    }
private:
    TaskData m_data;
};

#endif  // _TASKDATA_H_