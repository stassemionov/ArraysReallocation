#include <qralg.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

// Copy part of source matrix with double block layout
// to according part of destinating matrix with the same layout.
// Matrix part is specified by left upper corner coordinates.
void copy_minor_double_block(double* dst_ptr,
    const double* src_ptr,
    const int N,
    const int b1, const int b2,
    const int db1, const int db2,
    const int row_coord, const int col_coord)
{
    if (row_coord >= N || col_coord >= N)
    {
        return;
    }

    // Big blocks counts in matrix
    const int block_count_in_col = (int)ceil(1.0 * N / b1);
    const int block_count_in_row = (int)ceil(1.0 * N / b2);
    // Small blocks counts in big blocks
    const int dblock_count_in_bcol = (int)ceil(1.0 * b1 / db1);
    const int dblock_count_in_brow = (int)ceil(1.0 * b2 / db2);
    const int dblock_count_in_diff_bcol = (N % b1 == 0) ?
        dblock_count_in_bcol :
        (int)ceil(1.0 * (N % b1) / db1);
    const int dblock_count_in_diff_brow = (N % b2 == 0) ?
        dblock_count_in_brow :
        (int)ceil(1.0 * (N % b2) / db2);

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

// Fill part of source matrix with double block layout
// by specifed value 'val'.
// Matrix part is specified by left upper corner coordinates.
void set_minor_double_block(double* dst_ptr,
    const int val,
    const int N,
    const int b1, const int b2,
    const int db1, const int db2,
    const int row_coord, const int col_coord)
{
    if (row_coord >= N || col_coord >= N)
    {
        return;
    }

    // Big blocks counts in matrix
    const int block_count_in_col = (int) ceil(1.0 * N / b1);
    const int block_count_in_row = (int) ceil(1.0 * N / b2);
    // Small blocks counts in big blocks
    const int dblock_count_in_bcol = (int)ceil(1.0 * b1 / db1);
    const int dblock_count_in_brow = (int)ceil(1.0 * b2 / db2);
    const int dblock_count_in_diff_bcol = (N % b1 == 0) ?
        dblock_count_in_bcol :
        (int)ceil(1.0 * (N % b1) / db1);
    const int dblock_count_in_diff_brow = (N % b2 == 0) ?
        dblock_count_in_brow :
        (int)ceil(1.0 * (N % b2) / db2);

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

// Copy part of source matrix with block layout
// to according part of destinating matrix with the same layout.
// Matrix part is specified by left upper corner coordinates.
void copy_minor_block(double* dst_ptr,
    const double* src_ptr,
    const int N,
    const int b1, const int b2,
    const int row_coord, const int col_coord)
{
    if (row_coord >= N || col_coord >= N)
    {
        return;
    }

    // Big blocks counts in matrix
    const int block_count_in_col = (int)ceil(1.0 * N / b1);
    const int block_count_in_row = (int)ceil(1.0 * N / b2);
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
    const int copy_size_CORNER = (real_b2 - loc_col_coord) * sizeof(double);
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
        const int copy_size = jb_block_width * sizeof(double);
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
    const int copy_size_LEFT = (real_b2 - loc_col_coord) * sizeof(double);
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

// Fill part of source matrix with block layout by specifed value 'val'.
// Matrix part is specified by left upper corner coordinates.
void set_minor_block(double* dst_ptr,
    const int val,
    const int N,
    const int b1, const int b2,
    const int row_coord,
    const int col_coord)
{
    if (row_coord >= N || col_coord >= N)
    {
        return;
    }

    // Big blocks counts in matrix
    const int block_count_in_col = (int)ceil(1.0 * N / b1);
    const int block_count_in_row = (int)ceil(1.0 * N / b2);
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
    const int copy_size_CORNER = (real_b2 - loc_col_coord) * sizeof(double);
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
        const int copy_size = jb_block_width * sizeof(double);
        double* dst_row = dst_stripe + jb_shift * real_b1 +
            loc_row_coord * jb_block_width;
        for (int i = loc_row_coord; i < real_b1; ++i)
        {
            memset(dst_row, 0, copy_size);
            dst_row += jb_block_width;
        }
    }

    // * Copy of first big-block-column, which can be incomplete
    const int copy_size_LEFT = (real_b2 - loc_col_coord) * sizeof(double);
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

#ifdef __cplusplus
extern "C" {
#endif

double* QR_WY_tiled(double* A, const int N, const int b1, const int b2)
{    
    // Counts of blocks in matrix
    const int block_count_in_col = (int) ceil(1.0 * N / b1);
    const int block_count_in_row = (int) ceil(1.0 * N / b2);
    
    double* work_memory = (double*)(
            calloc(N + b2 + 2*N*b2 + 2*N*N, sizeof(double)));

    double* v =    work_memory;
    double* w =    v + N;
    double* W =    w + b2;
    double* Y =    W + N*b2;
    double* WY =   Y + N*b2;
    double* Abuf = WY + N*N;

    // Indices of current block column and stripe
    int curr_i_block_ind = -1;
    int curr_j_block_ind = -1;

    // 'lambda' is index of first column of current block-column
    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // Indices of big block, that contains element (lambda,lambda).
        // Remark: row-index value can increase with passing over diagonal
        //         when 'b1' != 'b2'.
        curr_i_block_ind = lambda / b1;
        ++curr_j_block_ind;
        // Vertical shift to diagonal block
        const int d1_shift = curr_i_block_ind * b1;
        // Index of first column of next block-column
        const int t = min(lambda + b2, N);
        // Real values of height of current block-stripe/column
        const int d2 = t - lambda;
        // Count of big blocks from the right of current big block in its block-stripe
        const int row_b_count = block_count_in_row - curr_j_block_ind - 1;
        // Count of big blocks below and including current big block in its block-column
        const int col_b_count = block_count_in_col - curr_i_block_ind;

        // \C2\FB\EF\EE\EB\ED\E5\ED\E8\E5 \EF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\E9 \ED\E0\E4 \F2\E5\EA\F3\F9\E8\EC \E1\EB\EE\EA\EE\EC
        for (int it = lambda; it < t; ++it)
        {
            const int loc_it = it - lambda;
            const int w_diff_val = t - it;
            const int A_diff_val = N - w_diff_val;
            // \C8\ED\E4\E5\EA\F1 \E1\EE\EB\FC\F8\EE\E3\EE \E1\EB\EE\EA\E0, \F1\EE\E4\E5\F0\E6\E0\F9\E5\E3\EE 'lambda+it'-\FB\E9 \E4\E8\E0\E3\EE\ED\E0\EB\FC\ED\FB\E9 \FD\EB\E5\EC\E5\ED\F2,
            // \F1\F7\E8\F2\E0\FF \EE\F2 \E1\EB\EE\EA\E0, \F1\EE\E4\E5\F0\E6\E0\F9\E5\E3\EE 'lambda'-\FB\E9 \FD\EB\E5\EC\E5\ED\F2
            const int b_begin = (it - d1_shift) / b1;

            memset(w, 0, d2*sizeof(double));

            double norm = 0.0;
            double scalar = 1.0;
            // * \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            const double* A_shifted_HV = A + it*N + it;
            double* v_shifted_HV = v + it;
            for (int j = 0; j < N - it; ++j)
            {
                const double buf = *A_shifted_HV;
                *v_shifted_HV++ = buf;
                norm += buf * buf;
                A_shifted_HV += N;
            }
            const double A_diag_el = A[it*N + it];
            const double diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + diag_el_sign * sqrt(norm);
            // \E8\F1\EF\EE\EB\FC\E7\F3\E5\F2\F1\FF \ED\EE\F0\EC\E8\F0\EE\E2\EA\E0, \F2.\F7. v[it] = 1
            for (int j = it + 1; j < N; ++j)
            {
                double buf = v[j] /= beta;
                scalar += buf * buf;
            }
            v[it] = 1.0;
            beta = -2.0 / scalar;

            // * \CF\F0\E8\EC\E5\ED\E5\ED\E8\E5 \EF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\FF

            // \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \E2\E5\EA\F2\EE\F0\E0 w = beta * (A(:,lambda:t-1)^t) * v .
            // \ED\E0\F7\E8\ED\E0\E5\EC \F1 it, \F2.\EA. it-\E0\FF \EC\E0\F2\F0\E8\F6\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0 \E4\E5\E9\F1\F2\E2\F3\E5\F2 \ED\E0 \EC\E0\F2\F0\E8\F6\F3 A(it:,it:)
            const int j_iter_count_TRANSF = t - it;
            for (int ib = b_begin; ib < col_b_count; ++ib)
            {
                const int ib_shift = d1_shift + ib*b1;
                const int lb = max(it, ib_shift);
                const int ub = min(ib_shift + b1, N);
                const double* A_row_shifted = A + lb*N + it;
                const double* v_shifted = v + lb;
                double* w_shifted = w + loc_it;

                for (int i = 0; i < ub - lb; ++i)
                {
                    const double vi = *v_shifted++;
                    for (int j = 0; j < j_iter_count_TRANSF; ++j)
                    {
                        (*w_shifted++) += (*A_row_shifted++) * vi;
                    }
                    A_row_shifted += A_diff_val;
                    w_shifted -= w_diff_val;
                }
            }
            for (int k = loc_it; k < d2; ++k)
            {
                w[k] *= beta;
            }

            // Transformation of current block-column:
            // A(it:N, loc_it:d2-1) += v[it:N] * (w[loc_it:d2-1]^t)
            for (int ib = b_begin; ib < col_b_count; ++ib)
            {
                const int mb_shift = d1_shift + ib * b1;
                const int i_lb = max(it, mb_shift);
                const int i_ub = min(mb_shift + b1, N);
                const double* v_shifted = v + i_lb;
                const double* w_shifted = w + loc_it;
                double* A_row_shifted = A + i_lb*N + it;

                for (int i = 0; i < i_ub - i_lb; ++i)
                {
                    const double vi = *v_shifted++;
                    for (int j = 0; j < j_iter_count_TRANSF; ++j)
                    {
                        (*A_row_shifted++) += vi * (*w_shifted++);
                    }
                    A_row_shifted += A_diff_val;
                    w_shifted -= w_diff_val;
                }
            }

            // * \C7\E0\EF\EE\EB\ED\E5\ED\E8\E5 \EF\EE\E4\E4\E8\E0\E3\ED\E0\EB\FC\ED\EE\E9 \F7\E0\F1\F2\E8 it-\E3\EE \F1\F2\EE\EB\E1\F6\E0 \EC\E0\F2\F0\E8\F6\FB A
            //   \E8\ED\F4\EE\F0\EC\E0\F2\E8\E2\ED\EE\E9 \F7\E0\F1\F2\FC\FE \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            double* A_shifted_fill = A + (it + 1)*N + it;
            const double* v_shifted_fill = v + it + 1;
            for (int j = it + 1; j < N; ++j)
            {
                (*A_shifted_fill) = (*v_shifted_fill++);
                A_shifted_fill += N;
            }
        }// FOR

         // * \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 WY-\EF\F0\E5\E4\F1\F2\E0\E2\EB\E5\ED\E8\FF \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF \EC\E0\F2\F0\E8\F6 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
        const int d = t - lambda;
        for (int it = 0; it < d; ++it)
        {
            const int shift = lambda + it;
            const int b_begin = (shift - d1_shift) / b1;

            memset(v + d1_shift, 0, (N - d1_shift)*sizeof(double));

            double scalar = 1, buf;
            // * \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \ED\EE\F0\EC\FB \E8 \F1\EA. \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            const double* A_shifted_HV = A + (shift + 1)*N + shift;
            double* v_shifted_HV = v + shift + 1;
            for (int j = shift + 1; j < N; ++j)
            {
                buf = *A_shifted_HV;
                scalar += buf * buf;
                (*v_shifted_HV++) = buf;
                A_shifted_HV += N;
            }
            v[shift] = 1.0;
            double beta = -2.0 / scalar;

            // it - \F7\E8\F1\EB\EE \ED\E0\EA\EE\EF\EB\E5\ED\ED\FB\F5 \F1\F2\EE\EB\E1\F6\EE\E2 \E2 \EC\E0\F2\F0\E8\F6\E0\F5 W \E8 Y (= \ED\EE\EC\E5\F0\F3 \E4\EE\E7\E0\EF\E8\F1\FB\E2\E0\E5\EC\EE\E3\EE \F1\F2\EE\EB\E1\F6\E0)
            if (it == 0)
            {
                // \CD\E0\F7\E0\EB\FC\ED\EE\E5 \E7\E0\EF\EE\EB\ED\E5\ED\E8\E5 \EF\E5\F0\E2\EE\E3\EE \F1\F2\EE\EB\E1\F6\E0
                // \E7\E0\EC\E5\F7\E0\ED\E8\E5: \EF\E5\F0\E2\FB\E5 lambda \F1\F2\F0\EE\EA \EC\E0\F2\F0\E8\F6 W \E8 Y - \ED\F3\EB\E5\E2\FB\E5,
                //            \EC\E0\F2\F0\E8\F6\FB \E8\EC\E5\FE\F2 \F1\F2\F3\EF\E5\ED\F7\E0\F2\FB\E9 \E2\E8\E4.
                double* Y_shifted_init = Y + lambda*b2;
                double* W_shifted_init = W + lambda*b2;
                const double* v_shifted_init = v + lambda;
                for (int i = 0; i < N - lambda; ++i)
                {
                    const double buf_v = *v_shifted_init++;
                    *Y_shifted_init = buf_v;
                    *W_shifted_init = beta * buf_v;
                    Y_shifted_init += b2;
                    W_shifted_init += b2;
                }
            }
            else
            {
                memset(w, 0, d * sizeof(double));

                const int W_Y_diff_val = b2 - it;

                // \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF (Y^t) * v
                double* wp_Yv = w;
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    const int i_lb = max(shift, d1_shift + ib * b1);
                    const int i_ub = min(d1_shift + (ib + 1) * b1, N);
                    const double* Y_row = Y + i_lb * b2;
                    const double* v_shifted = v + i_lb;

                    for (int i = 0; i < i_ub - i_lb; ++i)
                    {
                        const double vi = *v_shifted++;
                        for (int j = 0; j < it; ++j)
                        {
                            (*wp_Yv++) += (*Y_row++) * vi;
                        }
                        wp_Yv -= it;
                        Y_row += W_Y_diff_val;
                    }
                }

                // \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF W * ((Y^t) * v) <=> W * w
                const double* wp_Ww = w;
                for (int ib = 0; ib < col_b_count; ++ib)
                {
                    const int i_lb = d1_shift + ib * b1;
                    const int i_ub = min(i_lb + b1, N);
                    double* W_row = W + i_lb * b2;
                    double* Y_row_shifted = Y + i_lb * b2 + it;
                    double* v_shifted = v + i_lb;

                    for (int i = 0; i < i_ub - i_lb; ++i)
                    {
                        const double vi = *v_shifted++;
                        double sum = 0.0;
                        for (int j = 0; j < it; ++j)
                        {
                            sum += (*W_row++) * (*wp_Ww++);
                        }
                        *W_row = beta * (vi + sum);
                        *Y_row_shifted = vi;

                        wp_Ww -= it;
                        W_row += W_Y_diff_val;
                        Y_row_shifted += b2;
                    }
                }
            }
        }

        // * \CF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\E5 \EE\F1\F2\E0\EB\FC\ED\EE\E9 \F7\E0\F1\F2\E8 \EC\E0\F2\F0\E8\F6\FB A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        // Multiplication Y * W^t
#pragma omp parallel for schedule(dynamic, 1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib * b1;
            const int i_lb = max(lambda, i_mb_shift);
            const int i_ub = min(i_mb_shift + b1, N);
            const int i_iter_count = i_ub - i_lb;
            const int Y_diff_val = i_iter_count * b2;
            const int i_lb_mul_N = i_lb * N;
            const double* Y_row = Y + i_lb * b2;

            for (int jb = 0; jb < col_b_count; ++jb)
            {
                const int j_mb_shift = d1_shift + jb * b1;
                const int j_lb = max(lambda, j_mb_shift);
                const int j_ub = min(j_mb_shift + b1, N);
                const int j_iter_count = j_ub - j_lb;
                const int WY_diff_val = N - j_iter_count;
                const int W_diff_val = j_iter_count * b2;
                double* WY_row_shifted = WY + i_lb_mul_N + j_lb;
                const double* W_row = W + j_lb * b2;

                for (int i = 0; i < i_iter_count; ++i)
                {
                    const int kbound = min(d, i_lb + i - lambda + 1);
                    const int Wrow_diff_val = b2 - kbound;
                    for (int j = 0; j < j_iter_count; ++j)
                    {
                        double sum = 0.0;
                        for (int k = 0; k < kbound; ++k)
                        {
                            sum += (*Y_row++) * (*W_row++);
                        }
                        (*WY_row_shifted++) = sum;

                        Y_row -= kbound;
                        W_row += Wrow_diff_val;
                    }
                    Y_row += b2;
                    W_row -= W_diff_val;
                    WY_row_shifted += WY_diff_val;
                }
                Y_row -= Y_diff_val;
            }
        }

        const double* A_shifted_pre_copy = A + lambda * N + t;
        double* Abuf_shifted_pre_copy = Abuf + lambda * N + t;
        const int copy_length = (N - t) * sizeof(double);
        for (int i = lambda; i < N; ++i)
        {
            memcpy(Abuf_shifted_pre_copy, A_shifted_pre_copy, copy_length);
            Abuf_shifted_pre_copy += N;
            A_shifted_pre_copy += N;
        }

        // \D3\EC\ED\EE\E6\E5\ED\E8\E5 A(lamda:N-1, t:N-1) += (Y * W^t) * A(lamda:N-1, t:N-1)
#pragma omp parallel for schedule(dynamic, 1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_lb = lambda + ib * b1;
            const int i_ub = min(i_lb + b1, N);
            const int i_iter_count = i_ub - i_lb;
            const int Abuf_diff_val = i_iter_count * N;
            
            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_lb = t + jb * b2;
                const int j_ub = min(j_lb + b2, N);
                const int j_iter_count = j_ub - j_lb;
                const int Abuf_row_diff_val = N - j_iter_count;
                double* Abuf_row_shifted = Abuf + i_lb * N + j_lb;

                for (int kb = 0; kb < col_b_count; ++kb)
                {
                    const int k_lb = lambda + kb * b1;
                    const int k_ub = min(k_lb + b1, N);
                    const int k_iter_count = k_ub - k_lb;
                    const int A_row_diff_val = k_iter_count * N - 1;
                    const double* WY_row_shifted = WY + i_lb * N + k_lb;
                    const double* A_row_shifted = A + k_lb * N + j_lb;

                    for (int i = 0; i < i_iter_count; ++i)
                    {
                        for (int j = 0; j < j_iter_count; ++j)
                        {
                            double sum = *Abuf_row_shifted;
                            for (int k = 0; k < k_iter_count; ++k)
                            {
                                sum += (*WY_row_shifted++) * (*A_row_shifted);
                                A_row_shifted += N;
                            }
                            (*Abuf_row_shifted++) = sum;

                            WY_row_shifted -= k_iter_count;
                            A_row_shifted -= A_row_diff_val;
                        }
                        Abuf_row_shifted += Abuf_row_diff_val;
                        WY_row_shifted += N;
                        A_row_shifted -= j_iter_count;
                    }
                    Abuf_row_shifted -= Abuf_diff_val;
                }
            }
        }
                
        // \EA\EE\EF\E8\F0\EE\E2\E0\ED\E8\E5 \EF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\ED\EE\E9 \F7\E0\F1\F2\E8 \EC\E0\F2\F0\E8\F6\FB A
        double* A_shifted_post_copy = A + lambda * N + t;
        const double* Abuf_shifted_post_copy = Abuf + lambda * N + t;
        for (int i = lambda; i < N; ++i)
        {
            memcpy(A_shifted_post_copy, Abuf_shifted_post_copy, copy_length);
            Abuf_shifted_post_copy += N;
            A_shifted_post_copy += N;
        }

    }// WHILE

    free(work_memory);

    return A;
}

double* QR_WY_double_tiled(double* A,
    const int N, const int b1, const int b2, const int db1, const int db2)
{
    const int bsizes_ratio = (int)ceil(1.0 * b2 / b1);
    // Big blocks counts in matrix
    const int block_count_in_col = (int)ceil(1.0 * N / b1);
    const int block_count_in_row = (int)ceil(1.0 * N / b2);
    // Small blocks counts in big blocks
    const int dblock_count_in_bcol = (int)ceil(1.0 * b1 / db1);
    const int dblock_count_in_brow = (int)ceil(1.0 * b2 / db2);
    const int dblock_count_in_diff_bcol = (N % b1 == 0) ?
        dblock_count_in_bcol :
        (int)ceil(1.0 * (N % b1) / db1);
    const int dblock_count_in_diff_brow = (N % b2 == 0) ?
        dblock_count_in_brow :
        (int)ceil(1.0 * (N % b2) / db2);

    double* work_memory = (double*)(
            calloc(N + b2 + 2*N*b2 + 2*N*N + db1, sizeof(double)));

    double* v =       work_memory;
    double* w =       v + N;
    double* sum_vec = w + b2;
    double* W =       sum_vec + db1;
    double* Y =       W + N*b2;
    double* WY =      Y + N*b2;
    double* Abuf =    WY + N*N;

    int curr_i_block_ind = -1;
    int curr_j_block_ind = -1;

    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // Indices of big block, that contains element (lambda,lambda).
        // Remark: row-index value can increase with passing over diagonal
        //         when 'b1' != 'b2'.
        curr_i_block_ind = lambda / b1;
        ++curr_j_block_ind;

        // Index of first column of next block-column
        const int t = min(lambda + b2, N);
        // Count of big blocks from the right of current big block in its block-stripe
        const int row_b_count = block_count_in_row - curr_j_block_ind - 1;
        // Count of big blocks below and including current big block in its block-column
        const int col_b_count = block_count_in_col - curr_i_block_ind;
        // Horizontal and vertical shifts to diagonal block
        const int d1_shift = curr_i_block_ind * b1;
        const int d2_shift = curr_j_block_ind * b2;
        // Real values of height/width of current block-stripe/column
        const int d2 = t - lambda;
        // Count of small blocks in current big block in row-direction
        const int dblock_count_in_curr_block_row =
            (curr_j_block_ind == block_count_in_row - 1) ?
            dblock_count_in_diff_brow : dblock_count_in_brow;

        // * Execution of current block-column transformations...
        for (int it = 0; it < d2; ++it)
        {
            memset(w, 0, d2 * sizeof(double));
            memset(v + lambda, 0, it * sizeof(double));

            // Shift for 'lambda+it'-th element of matrix column
            const int full_shift = lambda + it;
            // \C8\ED\E4\E5\EA\F1 \E1\EE\EB\FC\F8\EE\E3\EE \E1\EB\EE\EA\E0, \F1\EE\E4\E5\F0\E6\E0\F9\E5\E3\EE 'lambda+it'-\FB\E9 \E4\E8\E0\E3\EE\ED\E0\EB\FC\ED\FB\E9 \FD\EB\E5\EC\E5\ED\F2,
            // \F1\F7\E8\F2\E0\FF \EE\F2 \E1\EB\EE\EA\E0, \F1\EE\E4\E5\F0\E6\E0\F9\E5\E3\EE 'lambda'-\FB\E9 \FD\EB\E5\EC\E5\ED\F2
            const int b_begin = (full_shift - d1_shift) / b1;
            // Index of small-block-column,
            // which contains 'lambda+it'-th column
            const int dblock_col_index = it / db2;

            // * Householder vector computation...
            double norm = 0.0;
            double scalar = 1.0;
            const double* A_shifted_HV = A + full_shift * N + full_shift;
            double* v_shifted_HV = v + full_shift;
            for (int i = 0; i < N - full_shift; ++i)
            {
                double buf = *A_shifted_HV;
                v_shifted_HV[i] = buf;
                norm += buf * buf;
                A_shifted_HV += N;
            }
            const double A_diag_el = A[full_shift*N + full_shift];
            const double A_diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + A_diag_el_sign * sqrt(norm);
            // \E8\F1\EF\EE\EB\FC\E7\F3\E5\F2\F1\FF \ED\EE\F0\EC\E8\F0\EE\E2\EA\E0, \F2.\F7. v[it] = 1
            double* v_shifted_NORM = v + full_shift + 1;
            for (int j = 0; j < N-full_shift-1; ++j)
            {
                double buf = v_shifted_NORM[j] /= beta;
                scalar += buf * buf;
            }
            v[full_shift] = 1.0;
            beta = -2.0 / scalar;

            // * Transformation applying...

            // Computation of vector w = beta * (A(:,lambda:t-1)^t) * v .
            // Starting with 'full_shift'-th row and column of matrix,
            // because 'full_shift'-th Householder matrix really operates
            // with A(full_shift:*, full_shift:t-1).
            // Pattern of data access is as follows:
            // 1. Fix big block (ib)
            // 2. Fix small block stripe in fixed big block (id)
            // 3. Fix small block in fixed stripe (jd)
            // 4. Fix row in fixed small block (i)
            // 5. Produce passing of fixed row (j)
            for (int ib = b_begin; ib < col_b_count; ++ib)
            {
                const int mb_shift = d1_shift + ib * b1;
                const int block_lb = max(full_shift, mb_shift) - mb_shift;
                const int block_ub = min(mb_shift + b1, N);
                const int id_lb = (block_lb == 0) ? 0 : block_lb / db1;
                const int id_ub = (ib == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;

                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int i_bd_shift = mb_shift + id * db1;
                    const int i_lb = max(full_shift, i_bd_shift) - i_bd_shift;
                    const int i_ub = min(db1, block_ub - i_bd_shift);
                    const int i_iter_count = i_ub - i_lb;
                    const double* A_dstripe = A + i_bd_shift * N;
                    const double* v_id_shifted = v + i_bd_shift + i_lb;

                    for (int jd = dblock_col_index; jd < dblock_count_in_curr_block_row; ++jd)
                    {
                        const int loc_drow_shift = jd * db2;
                        const int j_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int j_ub = min(d2 - loc_drow_shift, db2);
                        const int j_iter_count = j_ub - j_lb;
                        const int A_dif_val = N - j_iter_count;
                        double* w_shifted = w + loc_drow_shift + j_lb;
                        const double* A_row = (A_dstripe + i_lb*N) +
                            lambda + loc_drow_shift + j_lb;

                        for (int i = 0; i < i_iter_count; ++i)
                        {
                            const double vi = *v_id_shifted++;
                            for (int j = j_lb; j < j_ub; ++j)
                            {
                                (*w_shifted++) += (*A_row++) * vi;
                            }

                            w_shifted -= j_iter_count;
                            A_row += A_dif_val;
                        }
                        v_id_shifted -= i_iter_count;
                    }
                }
            }
            for (int k = it; k < d2; ++k)
            {
                w[k] *= beta;
            }

            // Transformation of current block-column:
            // A(fullshift:N, it:t-1) += v[fullshift:N] * (w[it:t-1]^t)
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                const int mb_shift = d1_shift + jb * b1;
                const int block_lb = max(full_shift, mb_shift) - mb_shift;
                const int block_ub = min(mb_shift + b1, N);
                const int jd_lb = (block_lb == 0) ? 0 : block_lb / db1;
                const int jd_ub = (jb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;

                for (int jd = jd_lb; jd < jd_ub; ++jd)
                {
                    const int bd_shift = mb_shift + jd * db1;
                    const int j_lb = max(full_shift, bd_shift) - bd_shift;
                    const int j_ub = min(db1, block_ub - bd_shift);
                    const int j_iter_count = j_ub - j_lb;
                    const int v_diff_val = j_ub - j_lb;
                    const int j_lb_mult_N = j_lb * N;
                    double* A_dstripe = A + bd_shift * N;
                    const double* v_jd_shifted = v + bd_shift + j_lb;

                    for (int kd = dblock_col_index; kd < dblock_count_in_curr_block_row; ++kd)
                    {
                        const int loc_drow_shift = kd * db2;
                        const int k_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int dblock_w = min(db2, d2 - loc_drow_shift);
                        const int k_iter_count = dblock_w - k_lb;
                        const int A_diff_val = N - k_iter_count;
                        const double* w_shifted = w + loc_drow_shift + k_lb;
                        double* A_line = (A_dstripe + j_lb_mult_N) +
                            lambda + loc_drow_shift + k_lb;

                        for (int j = 0; j < j_iter_count; ++j)
                        {
                            const double vj = (*v_jd_shifted++);
                            for (int k = 0; k < k_iter_count; ++k)
                            {
                                (*A_line++) += vj * (*w_shifted++);
                            }
                            A_line += A_diff_val;
                            w_shifted -= k_iter_count;
                        }

                        v_jd_shifted -= v_diff_val;
                    }
                }
            }

            // * \C7\E0\EF\EE\EB\ED\E5\ED\E8\E5 \EF\EE\E4\E4\E8\E0\E3\ED\E0\EB\FC\ED\EE\E9 \F7\E0\F1\F2\E8 it-\E3\EE \F1\F2\EE\EB\E1\F6\E0 \EC\E0\F2\F0\E8\F6\FB A
            //   \E8\ED\F4\EE\F0\EC\E0\F2\E8\E2\ED\EE\E9 \F7\E0\F1\F2\FC\FE \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            double* A_shifted_fill = A + (full_shift + 1)*N + full_shift;
            const double* v_shifted_fill = v + full_shift + 1;
            for (int j = full_shift + 1; j < N; ++j)
            {
                (*A_shifted_fill) = (*v_shifted_fill++);
                A_shifted_fill += N;
            }
        }// FOR

        // * Computation of WY-representation of Householder matrices production...

        const int set_len = min(bsizes_ratio*b1, N - d1_shift)*b2*sizeof(double);
        memset(W + d1_shift*b2, 0, set_len);
        memset(Y + d1_shift*b2, 0, set_len);
        for (int it = 0; it < d2; ++it)
        {
            memset(w, 0, d2*sizeof(double));
            memset(v + lambda, 0, it*sizeof(double));

            const int full_shift = lambda + it;
            const int b_begin = (full_shift - d1_shift) / b1;
            const int dblock_col_index = it / db2;

            // \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \ED\EE\F0\EC\FB \E8 \F1\EA. \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            double scalar = 1.0;
            double* A_shifted_HV = A + (full_shift + 1)*N + full_shift;
            double* v_shifted_HV = v + full_shift + 1;
            for (int j = full_shift + 1; j < N; ++j)
            {
                double buf = *A_shifted_HV;
                scalar += buf * buf;
                (*v_shifted_HV++) = buf;
                A_shifted_HV += N;
            }
            v[full_shift] = 1.0;
            double beta = -2.0 / scalar;

            // 'it' is count of accumulated columns in W and Y matrices
            // (equals index of columns, which is added at current iteration)
            if (it == 0)
            {
                // Initial first column filling.
                // Remark: first 'lambda' of rows in matrices W and Y are filled with zero,
                //           next 'd2' of rows consist lower triangular matrix.
                double* Y_shifted = Y + lambda*b2;
                double* W_shifted = W + lambda*b2;
                const double* v_shifted_WYinit = v + lambda;
                for (int i = lambda; i < N; ++i)
                {
                    double buf_WYinit = *v_shifted_WYinit++;
                    *Y_shifted = buf_WYinit;
                    *W_shifted = beta * buf_WYinit;
                    Y_shifted += b2;
                    W_shifted += b2;
                }
            }
            else
            {
                // Computation of production w = (Y^t) * v .
                // Pattern of data access is as follows:
                // 1. Fix big block (ib)
                // 2. Fix small block stripe in fixed big block (id)
                // 3. Fix small block in fixed stripe (jd)
                // 4. Fix row in fixed small block (i) and element of vector 'v'
                // 5. Produce passing of fixed row of Y (j) and vector 'w',
                //    using fixed element of vector 'v'
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    const int mb_shift = d1_shift + ib * b1;
                    const int block_lb = max(full_shift, mb_shift) - mb_shift;
                    const int block_ub = min(mb_shift + b1, N);
                    const int id_lb = (block_lb == 0) ? 0 : (block_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;

                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int bd_shift = mb_shift + id * db1;
                        const int i_lb = max(full_shift, bd_shift) - bd_shift;
                        const int i_ub = min(db1, block_ub - bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int v_diff_val = i_ub - i_lb;
                        const int i_lb_mult_b2 = i_lb * b2;
                        const double* Ydstripe = Y + bd_shift*b2;
                        const double* v_id_shifted = v + bd_shift + i_lb;

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it - loc_drow_shift, db2);
                            const int Y_diff_val = b2 - j_ub;
                            const double* Y_line = (Ydstripe + i_lb_mult_b2) + loc_drow_shift;
                            double* w_shifted = w + loc_drow_shift;

                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                const double vi = (*v_id_shifted++);
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    (*w_shifted++) += (*Y_line++) * vi;
                                }
                                Y_line += Y_diff_val;
                                w_shifted -= j_ub;
                            }

                            v_id_shifted -= v_diff_val;
                        }
                    }
                }

                // Computation of production W * ((Y^t) * v) <=> W * w .
                // Pattern of data access is as follows:
                // 1. Fix big block (ib)
                // 2. Fix small block stripe in fixed big block (id)
                // 3. Fix small block in fixed stripe (jd)
                // 4. Fix row in fixed small block (i)
                // 5. Produce passing of fixed row of W (j) and vector 'w'
                const int set_len_Ww = db1 * sizeof(double);
                for (int ib = 0; ib < col_b_count; ++ib)
                {
                    const int mb_shift = d1_shift + ib * b1;
                    const int iblock_lb = max(lambda, mb_shift) - mb_shift;
                    const int iblock_ub = min(mb_shift + b1, N);
                    const int id_lb = (iblock_lb == 0) ? 0 : (iblock_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;

                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int bd_shift = mb_shift + id * db1;
                        const int i_lb = max(lambda, bd_shift) - bd_shift;
                        const int i_ub = min(db1, iblock_ub - bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int i_lb_mult_b2 = i_lb * b2;
                        double* Wdstripe = W + bd_shift * b2;
                        double* Ydstripe = Y + (Wdstripe - W);
                        // Vector to contain results of production of W's row and vector w
                        memset(sum_vec, 0, set_len_Ww);

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it, loc_drow_shift + db2) - loc_drow_shift;
                            const int diff_val = b2 - j_ub;
                            const double* w_shifted = w + loc_drow_shift;
                            const double* Wrow = Wdstripe + i_lb_mult_b2 + loc_drow_shift;
                            double* sum_vec_shifted = sum_vec + i_lb;

                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                double sum_loc = 0.0;
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    sum_loc += (*Wrow++) * (*w_shifted++);
                                }
                                // Scalar production of W's row and vector w
                                (*sum_vec_shifted++) += sum_loc;

                                Wrow += diff_val;
                                w_shifted -= j_ub;
                            }
                        }
                        // Insertion of computed values into 'it'-th column of Y and W
                        const double* v_shifted_WYadd = v + bd_shift + i_lb;
                        // Pointers to 'it'-th column of Y and W
                        double* Wdblock_add = Wdstripe + i_lb_mult_b2 + it;
                        double* Ydblock_add = Ydstripe + (Wdblock_add - Wdstripe);
                        const double* sum_vec_shifted = sum_vec + i_lb;
                        for (int i = 0; i < i_iter_count; ++i)
                        {
                            double buf = v_shifted_WYadd[i];
                            *Wdblock_add = beta * (buf + sum_vec_shifted[i]);
                            *Ydblock_add = buf;
                            Wdblock_add += b2;
                            Ydblock_add += b2;
                        }
                    }
                }
            }
        }

        // * Transformation of remaining block-columns of A...
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        double* WY_shifted = WY + d1_shift*N + d1_shift;
        const int set_length = (N - d1_shift)*sizeof(double);
        for (int i = d1_shift; i < N; ++i)
        {
            memset(WY_shifted, 0, set_length);
            WY_shifted += N;
        }

        // Multiplication Y * W^t.
#pragma omp parallel for schedule(dynamic, 1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int iblock_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int iblock_ub = min(i_mb_shift + b1, N);
            const int id_lb = (iblock_lb == 0) ? 0 : (iblock_lb / db1);
            const int id_ub = (ib == col_b_count - 1) ?
                dblock_count_in_diff_bcol :
                dblock_count_in_bcol;

            for (int jb = 0; jb < col_b_count; ++jb)
            {
                const int j_mb_shift = d1_shift + jb*b1;
                const int jblock_lb = max(lambda, j_mb_shift) - j_mb_shift;
                const int jblock_ub = min(j_mb_shift + b1, N);
                const int jd_lb = (jblock_lb == 0) ? 0 : (jblock_lb / db1);
                const int jd_ub = (jb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;

                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int i_bd_shift = i_mb_shift + id * db1;
                    const int i_lb = max(lambda, i_bd_shift) - i_bd_shift;
                    const int i_ub = min(db1, iblock_ub - i_bd_shift);
                    const int i_iter_count = i_ub - i_lb;
                    const int i_lb_mul_b2 = i_lb * b2;
                    const int i_lb_mul_N = i_lb * N;
                    const double* Ydstripe = Y + i_bd_shift*b2;
                    double* WYdstripe = WY + i_bd_shift*N;

                    // Upper bound for 'kd'-level small-block-cycle is defined by
                    // lower triangular form of matrix Y.
                    // If 'id'-th small-block-stripe is crossing matrix diagonal,
                    // then upper bound for small-block-stripe passing equals
                    // index of small block, that is containing diagonal element
                    // with maximum index within 'id'-th small-block-stripe height bounds.
                    const int kd_ub = (i_bd_shift < t) ?
                        ((min(i_bd_shift + db1, iblock_ub) - lambda) / db2 + 1) :
                        dblock_count_in_curr_block_row;

                    for (int jd = jd_lb; jd < jd_ub; ++jd)
                    {
                        const int j_bd_shift = j_mb_shift + jd * db1;
                        const int j_lb = max(lambda, j_bd_shift) - j_bd_shift;
                        const int j_ub = min(db1, jblock_ub - j_bd_shift);
                        const int j_lb_mul_b2 = j_lb * b2;
                        const int WY_diff_val = N - (j_ub - j_lb);
                        const double* Wdstripe = W + j_bd_shift*b2;

                        for (int kd = 0; kd < kd_ub; ++kd)
                        {
                            const int k_loc_shift = kd * db2;
                            const int k_ub = min(d2 - k_loc_shift, db2);
                            const int W_diff_val = b2 - k_ub;

                            // Pointers to first row of small tile of W, Y and WY.
                            // Declared there because maust be updated on each iteration of 'kd'.
                            double* WYrow = (WYdstripe + i_lb_mul_N) + j_bd_shift + j_lb;
                            const double* Yrow = (Ydstripe + i_lb_mul_b2) + k_loc_shift;
                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                const double* Wrow = (Wdstripe + j_lb_mul_b2) + k_loc_shift;
                                for (int j = j_lb; j < j_ub; ++j)
                                {
                                    double sum = 0.0;
                                    for (int k = 0; k < k_ub; ++k)
                                    {
                                        sum += (*Yrow++) * (*Wrow++);
                                    }
                                    (*WYrow++) += sum;

                                    Yrow -= k_ub;
                                    Wrow += W_diff_val;
                                }
                                Yrow += b2;
                                WYrow += WY_diff_val;
                            }
                        }
                    }
                }
            }
        }

        const double* A_shifted = A + d1_shift*N + d2_shift;
        double* Abuf_shifted = Abuf + (A_shifted - A);
        const int copy_length = (N - d2_shift) * sizeof(double);
        for (int i = d1_shift; i < N; ++i)
        {
            memcpy(Abuf_shifted, A_shifted, copy_length);
            A_shifted += N;
            Abuf_shifted += N;
        }

        // Multiplication A(lamda:N-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1) <=>
        //                A(lamda:N-1, t:N-1) += WY * A(lamda:N-1, t:N-1) .
        // Layout pattern: (b1 x b2, d1 x d2) = (b1 x b1, d1 x d1) * (b1 x b2, d1 x d2).
        // Result of multiplication is added to temporal buffer 'Abuf',
        // which is already contain elements of 'A'.
        // It allows to avoid using of elements of matrix 'A',
        // which have results of multiplication, when we need their old values.
#pragma omp parallel for schedule(dynamic, 1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int iblock_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int iblock_ub = min(i_mb_shift + b1, N);
            const int id_lb = (iblock_lb == 0) ? 0 : (iblock_lb / db1);
            const int id_ub = (ib == col_b_count - 1) ?
                dblock_count_in_diff_bcol :
                dblock_count_in_bcol;

            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_mb_shift = t + jb*b2;
                const int jblock_ub = min(j_mb_shift + b2, N);
                const int jd_ub = (jb == row_b_count - 1) ?
                    dblock_count_in_diff_brow :
                    dblock_count_in_brow;

                for (int kb = 0; kb < col_b_count; ++kb)
                {
                    const int k_mb_shift = d1_shift + kb*b1;
                    const int kblock_lb = max(lambda, k_mb_shift) - k_mb_shift;
                    const int kblock_ub = min(k_mb_shift + b1, N);
                    const int kd_lb = (kblock_lb == 0) ? 0 : (kblock_lb / db1);
                    const int kd_ub = (kb == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;

                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int i_bd_shift = i_mb_shift + id * db1;
                        const int i_lb = max(lambda, i_bd_shift) - i_bd_shift;
                        const int i_ub = min(db1, iblock_ub - i_bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int i_lb_mul_N = i_lb * N;
                        double* Abuf_dstripe = Abuf + i_bd_shift*N;
                        const double* WY_dstripe = WY + (Abuf_dstripe - Abuf);

                        for (int jd = 0; jd < jd_ub; ++jd)
                        {
                            const int j_bd_shift = j_mb_shift + jd * db2;
                            const int j_ub = min(j_bd_shift + db2, jblock_ub) - j_bd_shift;
                            const int Abuf_dif_value = N - j_ub;
                            // Points to first row of ['id','jd'] small tile of 'Abuf'
                            double* Abuf_dstripe_shifted = (Abuf_dstripe + i_lb_mul_N) + j_bd_shift;

                            for (int kd = kd_lb; kd < kd_ub; ++kd)
                            {
                                const int k_bd_shift = k_mb_shift + kd * db1;
                                const int k_ub = min(k_bd_shift + db1, kblock_ub) - k_bd_shift;
                                const int A_diff_val = N * k_ub - 1;
                                // Points to 'i_lb'-th row of ['id','kd'] small tile of 'WY'
                                const double* WY_dstripe_shifted = (WY_dstripe + i_lb_mul_N) + k_bd_shift;
                                // Points to first row of ['kd','jd'] small tile of 'A'
                                const double* A_row = (A + k_bd_shift*N) + j_bd_shift;

                                double* Abuf_row = Abuf_dstripe_shifted;
                                for (int i = 0; i < i_iter_count; ++i)
                                {
                                    const double* WY_row = WY_dstripe_shifted;
                                    for (int j = 0; j < j_ub; ++j)
                                    {
                                        double sum = 0.0;
                                        for (int k = 0; k < k_ub; ++k)
                                        {
                                            sum += (*WY_row++) * (*A_row);
                                            A_row += N;
                                        }
                                        (*Abuf_row++) += sum;
                                        
                                        WY_row -= k_ub;
                                        A_row -= A_diff_val;
                                    }
                                    Abuf_row += Abuf_dif_value;
                                    WY_dstripe_shifted += N;
                                    A_row -= j_ub;
                                }
                            }
                        }
                    }
                }
            }
        }

        double* A_shifted_toUpdate = A + d1_shift*N + d2_shift;
        const double* Abuf_shifted_toUpdate = Abuf + (A_shifted_toUpdate - A);
        for (int i = d1_shift; i < N; ++i)
        {
            memcpy(A_shifted_toUpdate, Abuf_shifted_toUpdate, copy_length);
            A_shifted_toUpdate += N;
            Abuf_shifted_toUpdate += N;
        }
    }// WHILE

    free(work_memory);

    return A;
}

double* QR_WY_standard(double* A, const int N, const int r)
{
    double* work_memory = (double*)(
        calloc(N + r + 2*N*r + 2*N*N, sizeof(double)));

    double* v    = work_memory;
    double* w    = v + N;
    double* W    = w + r;
    double* Y    = W + N * r;
    double* WY   = Y + N * r;
    double* Abuf = WY + N * N;

    // lambda \F3\EA\E0\E7\FB\E2\E0\E5\F2 \ED\E0 \ED\E0\F7\E0\EB\EE \F2\E5\EA\F3\F9\E5\E3\EE \E1\EB\EE\EA\E0
    for (int lambda = 0; lambda < N; lambda += r)
    {
        // \D3\EA\E0\E7\FB\E2\E0\E5\F2 \ED\E0 \ED\E0\F7\E0\EB\EE \F1\EB\E5\E4\F3\FE\F9\E5\E3\EE \E1\EB\EE\EA\E0 \E7\E0 \F2\E5\EA\F3\F9\E8\EC
        const int t = min(lambda + r, N);

        // \C2\FB\EF\EE\EB\ED\E5\ED\E8\E5 \EF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\E9 \ED\E0\E4 \F2\E5\EA\F3\F9\E8\EC \E1\EB\EE\EA\EE\EC
        for (int it = lambda; it < t; ++it)
        {
            double norm = 0.0;
            double scalar = 1.0;
            double buf;
            // * \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            for (int j = it; j < N; ++j)
            {
                buf = A[j*N + it];
                v[j] = buf;
                norm += buf * buf;
            }
            const double A_diag_el = A[it * N + it];
            const double diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + diag_el_sign * sqrt(norm);
            // \E8\F1\EF\EE\EB\FC\E7\F3\E5\F2\F1\FF \ED\EE\F0\EC\E8\F0\EE\E2\EA\E0, \F2.\F7. v[it] = 1
            for (int j = it + 1; j < N; ++j)
            {
                buf = (v[j] /= beta);
                scalar += buf * buf;
            }
            v[it] = 1.0;

            // * \CF\F0\E8\EC\E5\ED\E5\ED\E8\E5 \EF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\FF
            beta = -2.0 / scalar;
            // \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \E2\E5\EA\F2\EE\F0\E0 w = beta * (A(:,lambda:t-1)^t) * v .
            // \ED\E0\F7\E8\ED\E0\E5\EC \F1 it, \F2.\EA. it-\E0\FF \EC\E0\F2\F0\E8\F6\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0 \E4\E5\E9\F1\F2\E2\F3\E5\F2 \ED\E0 \EC\E0\F2\F0\E8\F6\F3 A(it:,it:)
            for (int j = it; j < t; ++j)    
            {
                double sum = 0.0;
                for (int k = it; k < N; ++k)
                {
                    sum += A[k*N + j] * v[k];
                }
                w[j-it] = beta * sum;
            }
            // \CF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\E5 \F2\E5\EA\F3\F9\E5\E3\EE \E1\EB\EE\EA\E0 \EC\E0\F2\F0\E8\F6\FB A
            for (int j = it; j < N; ++j)
            {
                double* Aj = A + j*N;
                const double vj = v[j];
                for (int k = it; k < t; ++k)
                {
                    Aj[k] += vj * w[k-it];
                }
            }

            // * \C7\E0\EF\EE\EB\ED\E5\ED\E8\E5 \EF\EE\E4\E4\E8\E0\E3\ED\E0\EB\FC\ED\EE\E9 \F7\E0\F1\F2\E8 it-\E3\EE \F1\F2\EE\EB\E1\F6\E0 \EC\E0\F2\F0\E8\F6\FB A
            //   \E8\ED\F4\EE\F0\EC\E0\F2\E8\E2\ED\EE\E9 \F7\E0\F1\F2\FC\FE \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            for (int j = it + 1; j < N; ++j)
            {
                A[j*N + it] = v[j];
            }
        }// FOR
        
        // * \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 WY-\EF\F0\E5\E4\F1\F2\E0\E2\EB\E5\ED\E8\FF \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF \EC\E0\F2\F0\E8\F6 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
        const int d = t - lambda;
        memset(Y + lambda*r, 0, r * min(N - lambda, r) * sizeof(double));
        for (int it = 0; it < d; ++it)
        {
            const int shift = lambda + it;
            double scalar = 1.0, buf;
            memset(v + lambda, 0, it*sizeof(double));
            // * \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \ED\EE\F0\EC\FB \E8 \F1\EA. \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            for (int j = shift + 1; j < N; ++j)
            {
                buf = A[j*N + shift];
                scalar += buf * buf;
                v[j] = buf;
            }
            v[shift] = 1.0;
            double beta = -2.0 / scalar;

            // it - \F7\E8\F1\EB\EE \ED\E0\EA\EE\EF\EB\E5\ED\ED\FB\F5 \F1\F2\EE\EB\E1\F6\EE\E2 \E2 \EC\E0\F2\F0\E8\F6\E0\F5 W \E8 Y (= \ED\EE\EC\E5\F0\F3 \E4\EE\E7\E0\EF\E8\F1\FB\E2\E0\E5\EC\EE\E3\EE \F1\F2\EE\EB\E1\F6\E0)
            if (it == 0)
            {
                // \CD\E0\F7\E0\EB\FC\ED\EE\E5 \E7\E0\EF\EE\EB\ED\E5\ED\E8\E5 \EF\E5\F0\E2\EE\E3\EE \F1\F2\EE\EB\E1\F6\E0
                // \E7\E0\EC\E5\F7\E0\ED\E8\E5: \EF\E5\F0\E2\FB\E5 lambda \F1\F2\F0\EE\EA \EC\E0\F2\F0\E8\F6 W \E8 Y - \ED\F3\EB\E5\E2\FB\E5,
                //            \EC\E0\F2\F0\E8\F6\FB \E8\EC\E5\FE\F2 \F1\F2\F3\EF\E5\ED\F7\E0\F2\FB\E9 \E2\E8\E4.
                for (int i = lambda; i < N; ++i)
                {
                    Y[i*r] = v[i];
                    W[i*r] = beta * v[i];
                }
            }
            else
            {
                // \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF (Y^t) * v
                for (int i = 0; i < it; ++i)
                {
                    double sum = 0.0;
                    for (int j = shift; j < N; ++j)
                    {
                        sum += Y[j*r + i] * v[j];
                    }
                    w[i] = sum;
                }

                // \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF W * ((Y^t) * v)
                for (int i = lambda; i < N; ++i)
                {
                    double* Wi = W + i*r;
                    double sum = 0.0;
                    for (int j = 0; j < it; ++j)
                    {
                        sum += Wi[j] * w[j];
                    }

                    buf = v[i];
                    Wi[it] = beta * (buf + sum);
                    Y[Wi + it - W] = buf;
                }
            }
        }

        // * \CF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\E5 \EE\F1\F2\E0\EB\FC\ED\EE\E9 \F7\E0\F1\F2\E8 \EC\E0\F2\F0\E8\F6\FB A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)
        
        // \D3\EC\ED\EE\E6\E5\ED\E8\E5 Y * W^t
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = lambda; i < N; ++i)
        {
            double* WYi = WY + i*N;
            const double* Yi = Y + i*r;
            for (int j = lambda; j < N; ++j)
            {
                const double* Wj = W + j*r;
                double sum = 0.0;
                for (int k = 0; k < d; ++k)
                {
                    sum += Yi[k] * Wj[k];
                }
                WYi[j] = sum;
            }
        }

        // \D3\EC\ED\EE\E6\E5\ED\E8\E5 A(lamda:M-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1)
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = lambda; i < N; ++i)
        {
            double* Abuf_i = Abuf + i*N;
            const double* Ai = A + i*N;
            const double* WYi = WY + i*N;
            for (int j = t; j < N; ++j)
            {
                double sum = Ai[j];
                for (int k = lambda; k < N; ++k)
                {
                    sum += WYi[k] * A[k*N+j];
                }
                Abuf_i[j] = sum;
            }
        }

        // \EA\EE\EF\E8\F0\EE\E2\E0\ED\E8\E5 \EF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\ED\EE\E9 \F7\E0\F1\F2\E8 \EC\E0\F2\F0\E8\F6\FB A
        for (int i = lambda; i < N; ++i)
        {
            memcpy(A + i*N + t,
                   Abuf + i*N + t,
                   (N - t) * sizeof(double));
        }
    }// WHILE

    free(work_memory);

    return A;
}

double* QR_WY_block(double* A, const int N, const int b1, const int b2)
{
    const int block_count_in_col = (int)ceil(1.0 * N / b1);
    const int block_count_in_row = (int)ceil(1.0 * N / b2);

    double* work_memory = (double*)(
        calloc(N + b2 + 2*N*b2 + 2*N*N, sizeof(double)));

    double* v =    work_memory;
    double* w =    v + N;
    double* W =    w + b2;
    double* Y =    W + N*b2;
    double* WY =   Y + N*b2;
    double* Abuf = WY + N*N;

    int curr_i_block_ind = -1;
    int curr_j_block_ind = -1;
    // lambda \F3\EA\E0\E7\FB\E2\E0\E5\F2 \ED\E0 \ED\E0\F7\E0\EB\EE \F2\E5\EA\F3\F9\E5\E3\EE \E1\EB\EE\EA\E0
    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // \EA\EE\EE\F0\E4\E8\ED\E0\F2\FB \E1\EB\EE\EA\E0, \E2 \EA\EE\F2\EE\F0\EE\EC \ED\E0\F5\EE\E4\E8\F2\F1\FF \FD\EB\E5\EC\E5\ED\F2 (lambda,lambda).
        // \E7\E0\EC\E5\F7\E0\ED\E8\E5: i-\EA\EE\EE\F0\E4\E8\ED\E0\F2\E0 \EC\EE\E6\E5\F2 \E8\E7\EC\E5\ED\E8\F2\FC\F1\FF \E8\E7-\E7\E0 \F0\E0\E7\ED\FB\F5 \F0\E0\E7\EC\E5\F0\EE\E2 b1 \E8 b2 \E1\EB\EE\EA\EE\E2. 
        curr_i_block_ind = lambda / b1;
        ++curr_j_block_ind;

        // \D3\EA\E0\E7\FB\E2\E0\E5\F2 \ED\E0 \ED\E0\F7\E0\EB\EE \F1\EB\E5\E4\F3\FE\F9\E5\E3\EE \EF\EE \E3\EE\F0\E8\E7\EE\ED\F2\E0\EB\E8 \E1\EB\EE\EA\E0
        const int t = min(lambda + b2, N);
        // \EE\F1\F2\E0\E2\F8\E5\E5\F1\FF \F7\E8\F1\EB\EE \E1\EB\EE\EA\EE\E2 \F1\EF\F0\E0\E2\E0 \EE\F2 \F2\E5\EA\F3\F9\E5\E3\EE \EF\EE \F1\F2\F0\EE\EA\E5
        const int row_b_count = block_count_in_row - curr_j_block_ind - 1;
        // \EE\F1\F2\E0\E2\F8\E5\E5\F1\FF \F7\E8\F1\EB\EE \E1\EB\EE\EA\EE\E2 \E2\EA\EB\FE\F7\E0\FF \F2\E5\EA\F3\F9\E8\E9 \E8 \ED\E8\E6\E5 \EF\EE \F1\F2\EE\EB\E1\F6\F3
        const int col_b_count = block_count_in_col - curr_i_block_ind;
        // \F1\E4\E2\E8\E3\E8 \EF\EE \E2\E5\F0\F2\E8\EA\E0\EB\E8 \E8 \E3\EE\F0\E8\E7\EE\ED\F2\E0\EB\E8 \E4\EE \F2\E5\EA\F3\F9\E5\E3\EE \E1\EB\EE\EA\E0
        const int d1_shift = curr_i_block_ind * b1;
        const int d2_shift = curr_j_block_ind * b2;

        // \E4\E5\E9\F1\F2\E2\E8\F2\E5\EB\FC\ED\E0\FF \E2\E5\EB\E8\F7\E8\ED\E0 \F8\E8\F0\E8\ED\FB \F2\E5\EA\F3\F9\E5\E3\EE \E1\EB\EE\F7\ED\EE\E3\EE \F1\F2\EE\EB\E1\F6\E0
        const int d2 = t - lambda;
        
        // * \C2\FB\EF\EE\EB\ED\E5\ED\E8\E5 \EF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\E9 \ED\E0\E4 \F2\E5\EA\F3\F9\E8\EC \E1\EB\EE\F7\ED\FB\EC \F1\F2\EE\EB\E1\F6\EE\EC
        for (int it = 0; it < d2; ++it)
        {
            // \F1\E4\E2\E8\E3 \E4\EE 'lambda+it'-\E3\EE \FD\EB\E5\EC\ED\F2\E0 \E2 \F1\F2\EE\EB\E1\F6\E5
            const int full_shift = lambda + it;
            // \F4\E8\EA\F1\E0\F6\E8\FF \EF\E5\F0\E5\F5\EE\E4\E0 \E2 \F1\EB\E5\E4\F3\FE\F9\E8\E9 \EF\EE \E2\E5\F0\F2\E8\EA\E0\EB\E8 \E1\EB\EE\EA
            const int b_begin = (full_shift - d1_shift) / b1;

            memset(w, 0, d2 * sizeof(double));
            memset(v + lambda, 0, it * sizeof(double));
            
            // * \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            double norm = 0.0;
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                // \F1\E4\E2\E8\E3 \EB\E8\ED\E8\E8 \E1\EB\EE\EA\E0 \EE\F2 \EB\E8\ED\E8\E8 \EC\E0\F2\F0\E8\F6\FB
                const int mb_shift = d1_shift + jb * b1;
                // \E2\E5\F0\F5\ED\FF\FF \E8 \ED\E8\E6\ED\FF\FF \E3\F0\E0\ED\E8\F6\FB \F7\F2\E5\ED\E8\FF \F1\F2\F0\EE\EA \E2 \E1\EB\EE\EA\E5 jb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                const int j_lb = max(full_shift, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // \E2\FB\F1\EE\F2\E0 \E1\EB\EE\EA\E0 jb
                const int block_h = min(b1, N - mb_shift);
                // \F3\EA\E0\E7\E0\F2\E5\EB\FC \ED\E0 'it'-\FB\E9 \FD\EB\E5\EC\E5\ED\F2 'j_lb'-\EE\E9 \F1\F2\F0\EE\EA\E8 \E1\EB\EE\EA\E0 'jb'
                const double* Ablock_shifted = A + mb_shift*N +
                    block_h*d2_shift + j_lb*d2 + it;
                double* v_shifted = v + mb_shift + j_lb;
                // j - \ED\EE\EC\E5\F0 \F1\F2\F0\EE\EA\E8 \E2 \E1\EB\EE\EA\E5 jb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                for (int j = 0; j < j_ub - j_lb; ++j)
                {
                    double buf = *Ablock_shifted;
                    v_shifted[j] = buf;
                    norm += buf * buf;
                    Ablock_shifted += d2;
                }
            }

            const int diag_el_block_h = min(b1, N - (d1_shift + b_begin*b1));
            const double A_diag_el = *(A + d1_shift*N + b_begin*b1*N +
                d2_shift*diag_el_block_h + (full_shift % b1)*d2 + it);
            const double diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + diag_el_sign * sqrt(norm);
            // \E8\F1\EF\EE\EB\FC\E7\F3\E5\F2\F1\FF \ED\EE\F0\EC\E8\F0\EE\E2\EA\E0, \F2.\F7. v[it] = 1
            double scalar = 1.0;
            for (int j = full_shift + 1; j < N; ++j)
            {
                double buf = v[j] /= beta;
                scalar += buf * buf;
            }
            v[full_shift] = 1.0;
            beta = -2.0 / scalar;
        
            // * \CF\F0\E8\EC\E5\ED\E5\ED\E8\E5 \EF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\FF
            // \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \E2\E5\EA\F2\EE\F0\E0 w = beta * (A(:,lambda:t-1)^t) * v .
            // \CD\E0\F7\E8\ED\E0\E5\EC \F1 it, \F2.\EA. it-\E0\FF \EC\E0\F2\F0\E8\F6\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0 \E4\E5\E9\F1\F2\E2\F3\E5\F2 \ED\E0 \EC\E0\F2\F0\E8\F6\F3 A(it:,it:)
            for (int kb = b_begin; kb < col_b_count; ++kb)
            {
                // \F1\E4\E2\E8\E3 \EB\E8\ED\E8\E8 \E1\EB\EE\EA\E0 \EE\F2 \EB\E8\ED\E8\E8 \EC\E0\F2\F0\E8\F6\FB
                const int mb_shift = d1_shift + kb * b1;
                // \E2\E5\F0\F5\ED\FF\FF \E8 \ED\E8\E6\ED\FF\FF \E3\F0\E0\ED\E8\F6\FB \F7\F2\E5\ED\E8\FF \F1\F2\F0\EE\EA \E2 \E1\EB\EE\EA\E5 jb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                const int k_lb = max(full_shift, mb_shift) - mb_shift;
                const int k_ub = min(mb_shift + b1, N) - mb_shift;
                const int k_iter_count = k_ub - k_lb;
                // \C2\E5\EB\E8\F7\E8\ED\FB \F1\E4\E2\E8\E3\EE\E2 \F3\EA\E0\E7\E0\F2\E5\EB\E5\E9 \ED\E0 'A' \E8 \ED\E0 'v'
                // \F1 \E2\EE\E7\E2\F0\E0\F2\EE\EC \E8\F5 \ED\E0 \E8\F1\F5\EE\E4\ED\FB\E5 \EF\EE\E7\E8\F6\E8\E8
                // \E4\EB\FF \E2\FB\EF\EE\EB\ED\E5\ED\E8\FF \EE\F7\E5\F0\E5\ED\EE\E9 \E8\F2\E5\F0\E0\F6\E8\E8 \F3\EC\ED\EE\E6\E5\ED\E8\FF
                const int A_diff_val = k_iter_count * d2 - 1;
                // \E2\FB\F1\EE\F2\E0 \E1\EB\EE\EA\E0 jb
                const int block_h = min(b1, N - mb_shift);
                // \F3\EA\E0\E7\E0\F2\E5\EB\FC \ED\E0 \ED\E0\F7\E0\EB\EE \E1\EB\EE\EA\E0 jb
                const double* Ablock_shifted = A + mb_shift*N +
                    block_h*d2_shift + k_lb*d2 + it;
                const double* v_shifted = v + mb_shift + k_lb;
                // j - \ED\EE\EC\E5\F0 \F1\F2\EE\EB\E1\F6\E0 \E2 \E1\EB\EE\EA\E5 kb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                for (int j = it; j < d2; ++j)
                {
                    double sum = 0.0;
                    for (int k = 0; k < k_iter_count; ++k)
                    {
                        sum += (*Ablock_shifted) * v_shifted[k];
                        Ablock_shifted += d2;
                    }
                    w[j] += sum;

                    Ablock_shifted -= A_diff_val;
                }
            }
            for (int k = it; k < d2; ++k)
            {
                w[k] *= beta;
            }

            // \CF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\E5 \F2\E5\EA\F3\F9\E5\E3\EE \E1\EB\EE\EA\E0 \EC\E0\F2\F0\E8\F6\FB A
            const double* w_shifted_TRANSF = w + it;
            const int w_diff_val_TRANSF = d2 - it;
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                // \F1\E4\E2\E8\E3 \EB\E8\ED\E8\E8 \E1\EB\EE\EA\E0 \EE\F2 \EB\E8\ED\E8\E8 \EC\E0\F2\F0\E8\F6\FB
                const int mb_shift = d1_shift + jb*b1;
                // \E2\E5\F0\F5\ED\FF\FF \E8 \ED\E8\E6\ED\FF\FF \E3\F0\E0\ED\E8\F6\FB \F7\F2\E5\ED\E8\FF \F1\F2\F0\EE\EA \E2 \E1\EB\EE\EA\E5 jb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                const int j_lb = max(full_shift, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // \E2\FB\F1\EE\F2\E0 \E1\EB\EE\EA\E0 jb
                const int j_block_h = min(b1, N - mb_shift);
                // \F3\EA\E0\E7\E0\F2\E5\EB\FC \ED\E0 'it'-\FB\E9 \FD\EB\E5\EC\E5\ED\F2 'j_lb'-\EE\E9 \F1\F2\F0\EE\EA\E8 \E1\EB\EE\EA\E0 'jb'
                double* Arow_shifted = A + mb_shift*N +
                    j_block_h*d2_shift + j_lb*d2 + it;
                const double* v_shifted = v + mb_shift + j_lb;
                // j - \ED\EE\EC\E5\F0 \F1\F2\F0\EE\EA\E8 \E2 \E1\EB\EE\EA\E5 jb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                for (int j = 0; j < j_ub - j_lb; ++j)
                {
                    const double vj = v_shifted[j];
                    // k - \ED\EE\EC\E5\F0 \F1\F2\EE\EB\E1\F6\E0 \E2 \E1\EB\EE\EA\E5 jb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                    for (int k = 0; k < w_diff_val_TRANSF; ++k)
                    {
                        (*Arow_shifted++) += vj * (*w_shifted_TRANSF++);
                    }
                    Arow_shifted += it;
                    w_shifted_TRANSF -= w_diff_val_TRANSF;
                }
            }
            
            // * \C7\E0\EF\EE\EB\ED\E5\ED\E8\E5 \EF\EE\E4\E4\E8\E0\E3\ED\E0\EB\FC\ED\EE\E9 \F7\E0\F1\F2\E8 'it'-\E3\EE \F1\F2\EE\EB\E1\F6\E0 \EC\E0\F2\F0\E8\F6\FB 'A'
            //   \E8\ED\F4\EE\F0\EC\E0\F2\E8\E2\ED\EE\E9 \F7\E0\F1\F2\FC\FE \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            const int b_inc_begin = (full_shift+1 - curr_i_block_ind*b1) / b1;
            for (int jb = b_inc_begin; jb < col_b_count; ++jb)
            {
                // \F1\E4\E2\E8\E3 \EB\E8\ED\E8\E8 \E1\EB\EE\EA\E0 \EE\F2 \EB\E8\ED\E8\E8 \EC\E0\F2\F0\E8\F6\FB
                const int mb_shift = d1_shift + jb*b1;
                // \E2\E5\F0\F5\ED\FF\FF \E8 \ED\E8\E6\ED\FF\FF \E3\F0\E0\ED\E8\F6\FB \F7\F2\E5\ED\E8\FF \F1\F2\F0\EE\EA \E2 \E1\EB\EE\EA\E5 jb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                const int j_lb = max(full_shift+1, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // \E2\FB\F1\EE\F2\E0 \E1\EB\EE\EA\E0 jb
                const int j_block_h = min(b1, N - mb_shift);
                // \F3\EA\E0\E7\E0\F2\E5\EB\FC \ED\E0 \ED\E0\F7\E0\EB\EE \E1\EB\EE\EA\E0 jb
                double* Ablock_shifted = A + mb_shift*N +
                    j_block_h*d2_shift + j_lb*d2 + it;
                const double* v_shifted = v + mb_shift + j_lb;
                // j - \ED\EE\EC\E5\F0 \F1\F2\F0\EE\EA\E8 \E2 \E1\EB\EE\EA\E5 jb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                for (int j = 0; j < j_ub - j_lb; ++j)
                {
                    *Ablock_shifted = v_shifted[j];
                    Ablock_shifted += d2;
                }
            }
        }// FOR

        memset(W + lambda*b2, 0, (N - lambda)*b2*sizeof(double));
        memset(Y + lambda*b2, 0, (N - lambda)*b2*sizeof(double));

        // * \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 WY-\EF\F0\E5\E4\F1\F2\E0\E2\EB\E5\ED\E8\FF \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF \EC\E0\F2\F0\E8\F6 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
        for (int it = 0; it < d2; ++it)
        {
            // \F1\E4\E2\E8\E3 \E4\EE 'lambda+it'-\E3\EE \FD\EB\E5\EC\ED\F2\E0 \E2 \F1\F2\EE\EB\E1\F6\E5
            const int full_shift = lambda + it;
            // \F4\E8\EA\F1\E0\F6\E8\FF \EF\E5\F0\E5\F5\EE\E4\E0 \E2 \F1\EB\E5\E4\F3\FE\F9\E8\E9 \EF\EE \E2\E5\F0\F2\E8\EA\E0\EB\E8 \E1\EB\EE\EA
            const int b_begin = (full_shift - d1_shift) / b1;

            memset(w, 0, d2*sizeof(double));
            memset(v + lambda, 0, it*sizeof(double));

            double scalar = 1.0, buf = 0.0;
            // * \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \ED\EE\F0\EC\FB \E8 \F1\EA. \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF \E2\E5\EA\F2\EE\F0\E0 \D5\E0\F3\F1\F5\EE\EB\E4\E5\F0\E0
            for (int jb = (full_shift+1 - d1_shift)/b1; jb < col_b_count; ++jb)
            {
                // \F1\E4\E2\E8\E3 \EB\E8\ED\E8\E8 \E1\EB\EE\EA\E0 \EE\F2 \EB\E8\ED\E8\E8 \EC\E0\F2\F0\E8\F6\FB
                const int mb_shift = d1_shift + jb * b1;
                // \E2\E5\F0\F5\ED\FF\FF \E8 \ED\E8\E6\ED\FF\FF \E3\F0\E0\ED\E8\F6\FB \F7\F2\E5\ED\E8\FF \F1\F2\F0\EE\EA \E2 \E1\EB\EE\EA\E5 jb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                const int j_lb = max(full_shift+1, mb_shift) - mb_shift;
                const int j_ub = min(mb_shift + b1, N) - mb_shift;
                // \E2\FB\F1\EE\F2\E0 \E1\EB\EE\EA\E0 jb
                const int block_h = min(b1, N - mb_shift);
                // \F3\EA\E0\E7\E0\F2\E5\EB\FC \ED\E0 \ED\E0\F7\E0\EB\EE \E1\EB\EE\EA\E0 jb
                const double* Arow_shifted = A + mb_shift*N +
                    block_h*d2_shift + j_lb*d2 + it;
                double* v_shifted = v + mb_shift + j_lb;
                // j - \ED\EE\EC\E5\F0 \F1\F2\F0\EE\EA\E8 \E2 \E1\EB\EE\EA\E5 jb (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                for (int j = 0; j < j_ub - j_lb; ++j)
                {
                    buf = *Arow_shifted;
                    scalar += buf * buf;
                    v_shifted[j] = buf;
                    Arow_shifted += d2;
                }
            }
            v[full_shift] = 1.0;
            double beta = -2.0 / scalar;
            
            // 'it' - \F7\E8\F1\EB\EE \ED\E0\EA\EE\EF\EB\E5\ED\ED\FB\F5 \F1\F2\EE\EB\E1\F6\EE\E2 \E2 \EC\E0\F2\F0\E8\F6\E0\F5 'W' \E8 'Y'
            // (= \ED\EE\EC\E5\F0\F3 \E4\EE\EF\E8\F1\FB\E2\E0\E5\EC\EE\E3\EE \F1\F2\EE\EB\E1\F6\E0)
            if (it == 0)
            {
                // \CD\E0\F7\E0\EB\FC\ED\EE\E5 \E7\E0\EF\EE\EB\ED\E5\ED\E8\E5 \EF\E5\F0\E2\EE\E3\EE \F1\F2\EE\EB\E1\F6\E0
                // \E7\E0\EC\E5\F7\E0\ED\E8\E5: \EF\E5\F0\E2\FB\E5 'lambda' \F1\F2\F0\EE\EA \EC\E0\F2\F0\E8\F6 'W' \E8 'Y' - \ED\F3\EB\E5\E2\FB\E5,
                //            \EC\E0\F2\F0\E8\F6\E0 'Y' \E8\EC\E5\E5\F2 \F1\F2\F3\EF\E5\ED\F7\E0\F2\FB\E9 \E2\E8\E4.
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    // \F1\E4\E2\E8\E3 \EB\E8\ED\E8\E8 \E1\EB\EE\EA\E0 \EE\F2 \EB\E8\ED\E8\E8 \EC\E0\F2\F0\E8\F6\FB
                    const int mb_shift = d1_shift + ib * b1;
                    // \E2\E5\F0\F5\ED\FF\FF \E8 \ED\E8\E6\ED\FF\FF \E3\F0\E0\ED\E8\F6\FB \F7\F2\E5\ED\E8\FF \F1\F2\F0\EE\EA \E2 \E1\EB\EE\EA\E5 ib (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                    const int i_lb = max(full_shift, mb_shift) - mb_shift;
                    const int i_ub = min(mb_shift + b1, N) - mb_shift;
                    // \F3\EA\E0\E7\E0\F2\E5\EB\E8 \ED\E0 \ED\E0\F7\E0\EB\E0 \E1\EB\EE\EA\EE\E2 ib \EC\E0\F2\F0\E8\F6 Y \E8 W
                    double* Yrow = Y + (mb_shift + i_lb)*b2;
                    double* Wrow = W + (Yrow - Y);
                    const double* v_shifted = v + mb_shift + i_lb;
                    for (int i = 0; i < i_ub - i_lb; ++i)
                    {
                        double v_buf = v_shifted[i];
                        *Yrow = v_buf;
                        *Wrow = beta * v_buf;
                        Yrow += b2;
                        Wrow += b2;
                    }
                }
            }
            else
            {
                // \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF (Y^t) * v
                double* wp_Yv = w;
                const int W_Y_diff_val = b2 - it;
                for (int jb = b_begin; jb < col_b_count; ++jb)
                {
                    // \F1\E4\E2\E8\E3 \EB\E8\ED\E8\E8 \E1\EB\EE\EA\E0 \EE\F2 \EB\E8\ED\E8\E8 \EC\E0\F2\F0\E8\F6\FB
                    const int mb_shift = d1_shift + jb * b1;
                    // \E2\E5\F0\F5\ED\FF\FF \E8 \ED\E8\E6\ED\FF\FF \E3\F0\E0\ED\E8\F6\FB \F7\F2\E5\ED\E8\FF \F1\F2\F0\EE\EA \E2 \E1\EB\EE\EA\E5 jb (\E3\EB\EE\E1\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                    const int j_lb = max(full_shift, mb_shift);
                    const int j_ub = min(mb_shift + b1, N);
                    // \F3\EA\E0\E7\E0\F2\E5\EB\FC \ED\E0 \ED\E0\F7\E0\EB\EE \E1\EB\EE\EA\E0 jb
                    const double* Yrow = Y + j_lb*b2;
                    double* v_shifted = v + j_lb;
                    // j - \ED\EE\EC\E5\F0 \F1\F2\F0\EE\EA\E8 \E2 \E1\EB\EE\EA\E5 jb (\E3\EB\EE\E1\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                    for (int j = j_lb; j < j_ub; ++j)
                    {
                        const double vj = *v_shifted++;
                        for (int i = 0; i < it; ++i)
                        {
                            (*wp_Yv++) += (*Yrow++) * vj;
                        }
                        Yrow += W_Y_diff_val;
                        wp_Yv -= it;
                    }
                }
                
                // \C2\FB\F7\E8\F1\EB\E5\ED\E8\E5 \EF\F0\EE\E8\E7\E2\E5\E4\E5\ED\E8\FF W * ((Y^t) * v) <==> W * w
                const double* wp_Ww = w;
                for (int ib = 0; ib < col_b_count; ++ib)
                {
                    // \F1\E4\E2\E8\E3 \EB\E8\ED\E8\E8 \E1\EB\EE\EA\E0 \EE\F2 \EB\E8\ED\E8\E8 \EC\E0\F2\F0\E8\F6\FB
                    const int mb_shift = d1_shift + ib * b1;
                    // \E2\E5\F0\F5\ED\FF\FF \E8 \ED\E8\E6\ED\FF\FF \E3\F0\E0\ED\E8\F6\FB \F7\F2\E5\ED\E8\FF \F1\F2\F0\EE\EA \E2 \E1\EB\EE\EA\E5 ib (\EB\EE\EA\E0\EB\FC\ED\E0\FF \ED\F3\EC\E5\F0\E0\F6\E8\FF)
                    const int i_lb = max(lambda, mb_shift);
                    const int i_ub = min(mb_shift + b1, N);
                    // \F3\EA\E0\E7\E0\F2\E5\EB\E8 \ED\E0 \ED\E0\F7\E0\EB\E0 \E1\EB\EE\EA\EE\E2 ib \EC\E0\F2\F0\E8\F6 Y \E8 W
                    double* Yrow_shifted = Y + i_lb*b2 + it;
                    double* Wrow = W + i_lb*b2;
                    const double* v_shifted = v + i_lb;
                    for (int i = 0; i < i_ub - i_lb; ++i)
                    {
                        double sum = 0.0;
                        for (int j = 0; j < it; ++j)
                        {
                            sum += Wrow[j] * wp_Ww[j];
                        }
                        buf = v_shifted[i];
                        Wrow[it] = beta * (buf + sum);
                        *Yrow_shifted = buf;

                        Wrow += b2;
                        Yrow_shifted += b2;
                    }
                }
            }
        }

        // * \CF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\E8\E5 \EE\F1\F2\E0\EB\FC\ED\EE\E9 \F7\E0\F1\F2\E8 \EC\E0\F2\F0\E8\F6\FB A:
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        // \D3\EC\ED\EE\E6\E5\ED\E8\E5 Y * W^t
        // \CC\E0\F2\F0\E8\F6\E0 'WY' \E1\F3\E4\E5\F2 \E8\EC\E5\F2\FC 'b1' x 'b1' \F0\E0\E7\EC\E5\F9\E5\ED\E8\E5
#pragma omp parallel for schedule(dynamic, 1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib * b1;
            const int i_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int i_ub = min(b1, N - i_mb_shift);
            const int i_iter_count = i_ub - i_lb;
            const int Y_diff_val = (i_ub - i_lb) * b2;
            const int ib_stripe_h = min(b1, N - i_mb_shift);
            const double* Yrow = Y + (i_mb_shift + i_lb) * b2;
            double* WYstripe = WY + i_mb_shift*N + ib_stripe_h*d1_shift;

            for (int jb = 0; jb < col_b_count; ++jb)
            {
                const int j_mb_shift = d1_shift + jb * b1;
                const int j_lb = max(lambda, j_mb_shift) - j_mb_shift;
                const int j_ub = min(b1, N - j_mb_shift);
                const int j_iter_count = j_ub - j_lb;
                const int jb_block_w = min(b1, N - j_mb_shift);
                const int WY_diff_val = jb_block_w - (j_ub - j_lb);
                const int W_diff_val = (j_ub - j_lb) * b2;
                const double* Wrow = W + (j_mb_shift + j_lb) * b2;
                double* WYrow_shifted = WYstripe + jb * ib_stripe_h * b1 +
                    i_lb * jb_block_w + j_lb;

                for (int i = 0; i < i_iter_count; ++i)
                {
                    // \CE\E3\F0\E0\ED\E8\F7\E5\ED\E8\E5, \F1\E2\FF\E7\E0\ED\ED\EE\E5 \F1\EE \F1\F2\F3\EF\E5\ED\F7\E0\F2\FB\EC \E2\E8\E4\EE\EC \EC\E0\F2\F0\E8\F6\FB Y
                    const int kbound = min(t, i_mb_shift+i_lb + i+1) - lambda;
                    const int Wrow_diff_val = b2 - kbound;
                                        
                    for (int j = 0; j < j_iter_count; ++j)
                    {
                        double sum = 0.0;
                        for (int k = 0; k < kbound; ++k)
                        {
                            sum += (*Yrow++) * (*Wrow++);
                        }
                        (*WYrow_shifted++) = sum;

                        Yrow -= kbound;
                        Wrow += Wrow_diff_val;
                    }
                    WYrow_shifted += WY_diff_val;
                    Yrow += b2;
                    Wrow -= W_diff_val;
                }
                Yrow -= Y_diff_val;
            }
        }

        copy_minor_block(Abuf, A, N, b1, b2, d1_shift, d2_shift);

        // \D3\EC\ED\EE\E6\E5\ED\E8\E5 A(lamda:M-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1)
        // \D0\E0\E7\EC\E5\F9\E5\ED\E8\E5: (b1 x b2) = (b1 x b1) * (b1 x b2)
#pragma omp parallel for schedule(dynamic, 1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib * b1;
            const int i_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int ib_stripe_h = min(b1, N - i_mb_shift);
            const int i_iter_count = ib_stripe_h - i_lb;
            const double* A_stripe_shifted = A + i_mb_shift*N + t*ib_stripe_h;
            double* Abuf_stripe_shifted = Abuf + (A_stripe_shifted - A);
            const double* WY_stripe_shifted = WY + i_mb_shift*N + d1_shift*ib_stripe_h;

            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_mb_shift = t + jb * b2;
                const int jb_block_w = min(b2, N - j_mb_shift);
                const int Abuf_diff_val = (ib_stripe_h - i_lb) * jb_block_w;
                double* Abuf_row = Abuf_stripe_shifted +
                    jb * ib_stripe_h * b2 + i_lb * jb_block_w;

                for (int kb = 0; kb < col_b_count; ++kb)
                {
                    const int k_mb_shift = d1_shift + kb * b1;
                    const int k_lb = max(lambda, k_mb_shift) - k_mb_shift;
                    const int kb_block_size = min(b1, N - k_mb_shift);
                    const int k_iter_count = kb_block_size - k_lb;
                    const int Arow_diff_val = k_iter_count * jb_block_w - 1;
                    const double* WYrow_shifted = WY_stripe_shifted +
                        kb * ib_stripe_h * b1 + i_lb * kb_block_size + k_lb;
                    const double* Arow_shifted = A + k_mb_shift * N +
                        j_mb_shift * kb_block_size + k_lb * jb_block_w;

                    for (int i = 0; i < i_iter_count; ++i)
                    {
                        for (int j = 0; j < jb_block_w; ++j)
                        {
                            double sum = 0.0;
                            for (int k = 0; k < k_iter_count; ++k)
                            {
                                sum += WYrow_shifted[k] * (*Arow_shifted);
                                Arow_shifted += jb_block_w;
                            }
                            (*Abuf_row++) += sum;

                            Arow_shifted -= Arow_diff_val;
                        }
                        WYrow_shifted += kb_block_size;
                        Arow_shifted -= jb_block_w;
                    }
                    Abuf_row -= Abuf_diff_val;
                }
            }
        }

        // \CA\EE\EF\E8\F0\EE\E2\E0\ED\E8\E5 \EF\F0\E5\EE\E1\F0\E0\E7\EE\E2\E0\ED\ED\EE\E9 \F7\E0\F1\F2\E8 \EC\E0\F2\F0\E8\F6\FB A
        copy_minor_block(A, Abuf, N, b1, b2, d1_shift, d2_shift);

    }// WHILE

    free(work_memory);

    return A;
}

double* QR_WY_double_block(double* A,
    const int N, const int b1, const int b2, const int db1, const int db2)
{
    const int bsizes_ratio = (int)ceil(1.0 * b2 / b1);
    // Big blocks counts in matrix
    const int block_count_in_col = (int)ceil(1.0 * N / b1);
    const int block_count_in_row = (int)ceil(1.0 * N / b2);
    // Small blocks counts in big blocks
    const int dblock_count_in_bcol = (int)ceil(1.0 * b1 / db1);
    const int dblock_count_in_brow = (int)ceil(1.0 * b2 / db2);
    const int dblock_count_in_diff_bcol = (N % b1 == 0) ?
        dblock_count_in_bcol :
        (int)ceil(1.0 * (N % b1) / db1);
    const int dblock_count_in_diff_brow = (N % b2 == 0) ?
        dblock_count_in_brow :
        (int)ceil(1.0 * (N % b2) / db2);

    double* work_memory = (double*)(
            calloc(N + b2 + 2*N*b2 + 2*N*N + db1, sizeof(double)));

    double* v =       work_memory;
    double* w =       v + N;
    double* sum_vec = w + b2;
    double* W =       sum_vec + db1;
    double* Y =       W + N*b2;
    double* WY =      Y + N*b2;
    double* Abuf =    WY + N*N;

    int curr_i_block_ind = -1;
    int curr_j_block_ind = -1;
    
    for (int lambda = 0; lambda < N; lambda += b2)
    {
        // Indices of big block, that contains element (lambda,lambda).
        // Remark: row-index can be not constant when b1 != b2.
        curr_i_block_ind = lambda / b1;
        ++curr_j_block_ind;

        // Index of first column of next block-column
        const int t = min(lambda + b2, N);
        // Count of big blocks from the right of current big block in its block-stripe
        const int row_b_count = block_count_in_row - curr_j_block_ind - 1;
        // Count of big blocks below and including current big block in its block-column
        const int col_b_count = block_count_in_col - curr_i_block_ind;
        // Horizontal and vertical shifts to diagonal block
        const int d1_shift = curr_i_block_ind * b1;
        const int d2_shift = curr_j_block_ind * b2;
        // Real value of width of current block-column
        const int d2 = t - lambda;
        // Count of small blocks in current big block in row-direction
        const int dblock_count_in_curr_block_row = 
            (curr_j_block_ind == block_count_in_row - 1) ?
            dblock_count_in_diff_brow : dblock_count_in_brow;
        
        // * Execution of current block-column transformations...
        for (int it = 0; it < d2; ++it)
        {
            memset(w, 0, d2 * sizeof(double));
            memset(v + lambda, 0, it * sizeof(double));

            // Shift for 'lambda+it'-th element of matrix column
            const int full_shift = lambda + it;
            // \C8\ED\E4\E5\EA\F1 \E1\EE\EB\FC\F8\EE\E3\EE \E1\EB\EE\EA\E0, \F1\EE\E4\E5\F0\E6\E0\F9\E5\E3\EE 'lambda+it'-\FB\E9 \E4\E8\E0\E3\EE\ED\E0\EB\FC\ED\FB\E9 \FD\EB\E5\EC\E5\ED\F2,
            // \F1\F7\E8\F2\E0\FF \EE\F2 \E1\EB\EE\EA\E0, \F1\EE\E4\E5\F0\E6\E0\F9\E5\E3\EE 'lambda'-\FB\E9 \FD\EB\E5\EC\E5\ED\F2
            const int b_begin = (full_shift - d1_shift) / b1;
            // Index of small-block-column,
            // which contains 'lambda+it'-th column
            const int dblock_col_index = it / db2;
            // Index of small-block-row of diagonal big block,
            // which contains 'full_shift'-th row of matrix
            const int dblock_row_index = (full_shift - (d1_shift + b_begin*b1)) / db1;
            // Width of this small block
            const int dblock_width = min(db2, d2 - dblock_col_index*db2);
            const int loc_it = it - dblock_col_index*db2;

            // * Householder vector computation...
            // Using of 'it'-th column of current block-column,
            // starting with (lambda+it)-th position
            double norm = 0.0;
            for (int ib = b_begin; ib < col_b_count; ++ib)
            {
                // Shift of current big block first line
                // from matrix first line
                const int mb_shift = d1_shift + ib * b1;
                // Upper and lower bounds for 'ib'-th block lines reading.
                // 'i_lb' is required, because we need to use block rows
                // starting with 'full_shift'-th matrix row,
                // that can be mapped not only to first block row.
                const int i_lb = max(full_shift, mb_shift) - mb_shift;  // local indexation
                // 'i_ub' is required, because last matrix row,
                // that can be mapped not only to last block row, but smaller.
                const int i_ub = min(mb_shift + b1, N);  // global indexation
                // 'ib'-th block height (number of rows)
                const int ib_block_h = min(b1, N - mb_shift);
                // 'ib'-th block beginning pointer
                const double* Ablock = A + mb_shift * N + ib_block_h * d2_shift;

                // Upper and lower bounds for small blocks stripes reading.
                // Reasons of using these bounds are the same as 'i_lb' and 'i_ub'.
                const int id_lb = (i_lb == 0) ? 0 : (i_lb / db1);
                const int id_ub = (ib == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int loc_shift = id * db1;
                    // Shift of current small-block-row from matrix first row
                    const int bd_shift = mb_shift + loc_shift;
                    // Bounds for current small block reading.
                    // Reasons of using are the same as 'i_lb' and 'i_ub'
                    const int lb = max(full_shift, bd_shift) - bd_shift;
                    const int ub = min(bd_shift + db1, i_ub) - bd_shift;
                    // 'id'-th small block height (number of rows)
                    const int id_block_h = min(db1, ib_block_h - loc_shift);
                    // Pointer to beginning of required small block in 'id'-th stripe
                    const double* Adblock_shifted = Ablock + loc_shift*d2 +
                        dblock_col_index*id_block_h*db2 + lb*dblock_width + loc_it;
                    double* v_shifted = v + bd_shift + lb;
                    
                    for (int i = 0; i < ub - lb; ++i)
                    {
                        double buf = *Adblock_shifted;
                        v_shifted[i] = buf;
                        norm += buf * buf;

                        Adblock_shifted += dblock_width;
                    }
                }
            }
                        
            // Diagonal big block height
            const int diag_el_block_h = min(b1, N - (d1_shift + b_begin*b1));
            // Diagonal small block height
            const int diag_el_dblock_h = min(db1, diag_el_block_h - dblock_row_index*db1);
            const double A_diag_el = *(A + d1_shift*N + b_begin*b1*N +  // shift for big block stripe
                d2_shift * diag_el_block_h +  // shift for big block beginning in its stripe
                ((full_shift % b1) / db1)*db1*d2 +  // shift for small block stripe
                (dblock_col_index*diag_el_dblock_h*db2) +   // shift for small block beginning in its stripe
                ((full_shift % b1) % db1)*dblock_width + loc_it);  // shift for diagonal element in small block

            const double diag_el_sign = (A_diag_el < 0) ? -1 : 1;
            double beta = A_diag_el + diag_el_sign * sqrt(norm);
            // Use normalization such that v[it] = 1
            double scalar = 1.0;
            double* v_shifted_NORM = v + full_shift + 1;
            for (int j = 0; j < N - full_shift - 1; ++j)
            {
                double buf = v_shifted_NORM[j] /= beta;
                scalar += buf * buf;
            }
            v[full_shift] = 1.0;
            beta = -2.0 / scalar;

            // * Transformation applying...

            // Computation of vector w = beta * (A(:,lambda:t-1)^t) * v .
            // Starting with 'full_shift'-th row and column of matrix,
            // because 'full_shift'-th Householder matrix really operates
            // with A(full_shift:*, full_shift:t-1).
            // Pattern of data access is as follows:
            // 1. Fix big block (kb)
            // 2. Fix small block stripe in fixed big block (kd)
            // 3. Fix small block in fixed stripe (jd)
            // 4. Fix column in fixed small block (j)
            // 5. Produce passing of fixed column (k)
            for (int kb = b_begin; kb < col_b_count; ++kb)
            {
                const int mb_shift = kb*b1 + d1_shift;
                const int k_lb = max(full_shift, mb_shift) - mb_shift;
                const int k_ub = min(mb_shift + b1, N);
                const int block_h = min(b1, N - mb_shift);
                const double* Ablock = A + mb_shift*N + block_h*d2_shift;
                const double* v_kb_shift = v + mb_shift;
                
                const int kd_lb = (k_lb == 0) ? 0 : k_lb / db1;
                const int kd_ub = (kb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int kd = kd_lb; kd < kd_ub; ++kd)
                {
                    const int loc_shift = kd * db1;
                    const int bd_shift = loc_shift + mb_shift;
                    const int lb = max(full_shift, bd_shift) - bd_shift;
                    const int ub = min(db1, k_ub - bd_shift);
                    const int k_iter_count = ub - lb;
                    const int dblock_h = min(db1, block_h - loc_shift);
                    // Pointer to beginning of 'kb'-th small block stripe
                    const double* Adstripe = Ablock + loc_shift*d2;
                    const double* v_kd_shifted = v_kb_shift + loc_shift + lb;
                    
                    // const int& jd_lb = dblock_col_index;
                    // const int& jd_ub = dblock_count_in_curr_block_row;
                    for (int jd = dblock_col_index; jd < dblock_count_in_curr_block_row; ++jd)
                    {
                        const int loc_drow_shift = jd * db2;
                        // Width of this small block
                        const int dblock_w = min(db2, d2 - loc_drow_shift);
                        const int j_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int j_ub = min(d2 - loc_drow_shift, db2);
                        const int j_iter_count = j_ub - j_lb;
                        const int A_diff_val = k_iter_count * dblock_w - 1;
                        // Pointer to beginning of 'jb'-th small block
                        const double* Arow_shifted = Adstripe +
                            dblock_h * loc_drow_shift + lb * dblock_w + j_lb;
                        double* w_shifted = w + loc_drow_shift + j_lb;

                        for (int j = 0; j < j_iter_count; ++j)
                        {
                            double sum = 0.0;
                            for (int k = 0; k < k_iter_count; ++k)
                            {
                                sum += (*Arow_shifted) * v_kd_shifted[k];
                                Arow_shifted += dblock_w;
                            }
                            Arow_shifted -= A_diff_val;
                            w_shifted[j] += sum;
                        }
                    }
                }
            }
            for (int k = it; k < d2; ++k)
            {
                w[k] *= beta;
            }

            // Transformation of current block column of matrix A:
            // A(fullshift:N, it:t-1) += v[fullshift:N] * (w[it:t-1]^t)
            for (int jb = b_begin; jb < col_b_count; ++jb)
            {
                const int mb_shift = d1_shift + jb * b1;
                const int jb_lb = max(full_shift, mb_shift) - mb_shift;
                const int jb_ub = min(mb_shift + b1, N);
                const int jb_block_h = min(b1, N - mb_shift);
                double* Ablock = A + mb_shift*N + jb_block_h*d2_shift;

                const int jd_lb = (jb_lb == 0) ? 0 : jb_lb / db1;
                const int jd_ub = (jb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int jd = jd_lb; jd < jd_ub; ++jd)
                {
                    const int loc_shift = jd * db1;
                    const int bd_shift = loc_shift + mb_shift;
                    const int j_lb = max(full_shift, bd_shift) - bd_shift;
                    const int j_ub = min(db1, jb_ub - bd_shift);
                    const int j_it_count = j_ub - j_lb;
                    const int jd_block_h = min(db1, jb_block_h - loc_shift);
                    double* Adstripe = Ablock + loc_shift*d2;
                    const double* v_shifted = v + bd_shift + j_lb;

                    for (int kd = dblock_col_index; kd < dblock_count_in_curr_block_row; ++kd)
                    {
                        const int loc_drow_shift = kd * db2;
                        const int k_lb = max(it, loc_drow_shift) - loc_drow_shift;
                        const int dblock_w = min(db2, d2 - loc_drow_shift);
                        const int k_iter_count = dblock_w - k_lb;
                        double* Arow_shifted = Adstripe + jd_block_h * loc_drow_shift +
                            j_lb * dblock_w + k_lb;
                        const double* w_shifted = w + loc_drow_shift + k_lb;

                        for (int j = 0; j < j_it_count; ++j)
                        {
                            const double vj = v_shifted[j];
                            for (int k = 0; k < k_iter_count; ++k)
                            {
                                (*Arow_shifted++) += vj * (*w_shifted++);
                            }

                            w_shifted -= k_iter_count;
                            Arow_shifted += k_lb;
                        }
                    }
                }
            }

            // Filling up of underdiagonal part of matrix 'full_shift'-th column
            // with informative part of Householder vector 'v'
            const int b_inc_begin = (full_shift + 1 - d1_shift) / b1;
            for (int ib = b_inc_begin; ib < col_b_count; ++ib)
            {
                const int mb_shift = d1_shift + ib*b1;
                const int ib_lb = max(full_shift + 1, mb_shift) - mb_shift;
                const int ib_ub = min(mb_shift + b1, N);
                const int ib_block_h = min(b1, N - mb_shift);
                double* Ablock = A + mb_shift*N + ib_block_h*d2_shift;

                const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
                const int id_ub = (ib == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int loc_shift = id * db1;
                    const int bd_shift = loc_shift + mb_shift;
                    const int i_lb = max(full_shift + 1, bd_shift) - bd_shift;
                    const int i_ub = min(db1, ib_ub - bd_shift);
                    const int id_block_h = min(db1, ib_block_h - loc_shift);
                    double* Arow_shifted = Ablock + loc_shift*d2 +
                        dblock_col_index*id_block_h*db2 +
                        i_lb*dblock_width + loc_it;
                    const double* v_shifted = v + bd_shift + i_lb;
                    
                    for (int i = 0; i < i_ub - i_lb; ++i)
                    {
                        *Arow_shifted = v_shifted[i];
                        Arow_shifted += dblock_width;
                    }
                }
            }
        }// FOR

        // * Computation of WY-representation of Householder matrices production...
        memset(W + d1_shift*b2, 0, min(bsizes_ratio*b1, N - d1_shift)*b2*sizeof(double));
        memset(Y + d1_shift*b2, 0, min(bsizes_ratio*b1, N - d1_shift)*b2*sizeof(double));
        for (int it = 0; it < d2; ++it)
        {
            memset(w, 0, d2 * sizeof(double));
            memset(v + lambda, 0, it * sizeof(double));

            const int full_shift = lambda + it;
            const int b_begin = (full_shift - d1_shift) / b1;
            const int dblock_col_index = it / db2;
            const int dblock_width = min(db2, d2 - dblock_col_index * db2);
            // Separated value is required because W and Y matrices are
            // always storing as N x b2 with d1 x d2 small blocks. 
            // So, if actual width of column is less then 'b2',
            // then form of small blocks won't change.
            const int dblock_width_WY = min(db2, b2 - dblock_col_index*db2);
            const int loc_it = it - dblock_col_index*db2;

            double scalar = 1.0;
            const int b_inc_begin = (full_shift - d1_shift) / b1;
            for (int ib = b_inc_begin; ib < col_b_count; ++ib)
            {
                const int mb_shift = d1_shift + ib * b1;
                const int i_lb = max(full_shift + 1, mb_shift) - mb_shift;
                const int i_ub = min(mb_shift + b1, N);
                const int ib_block_h = min(b1, N - mb_shift);
                const double* Ablock = A + mb_shift * N +
                    ib_block_h * d2_shift;

                const int id_lb = (i_lb == 0) ? 0 : (i_lb / db1);
                const int id_ub = (ib == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int loc_shift = id * db1;
                    const int bd_shift = mb_shift + loc_shift;
                    const int lb = max(full_shift + 1, bd_shift) - bd_shift;
                    const int ub = min(db1, i_ub - bd_shift);
                    const int id_block_h = min(db1, ib_block_h - loc_shift);
                    const double* Adblock_shifted = Ablock + loc_shift*d2 +
                        dblock_col_index*id_block_h*db2 +
                        loc_it + lb*dblock_width;
                    double* v_shifted = v + bd_shift + lb;

                    for (int i = 0; i < ub - lb; ++i)
                    {
                        double buf = *Adblock_shifted;
                        Adblock_shifted += dblock_width;
                        scalar += buf * buf;
                        v_shifted[i] = buf;
                    }
                }
            }
            v[full_shift] = 1.0;
            double beta = -2.0 / scalar;
                        
            // 'it' is count of accumulated columns in W and Y matrices
            // (equals index of columns, which is added at current iteration)
            if (it == 0)
            {
                // Initial first column filling
                // Remark 1: first 'lambda' of rows in matrices W and Y are filled with zero,
                //           next 'd2' of rows consist lower triangular matrix.
                // Remark 2: W and Y will be allocated with (B1,B2,D1,D2)-layout.
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    const int mb_shift = ib*b1 + d1_shift;
                    const int i_lb = max(lambda, mb_shift) - mb_shift;
                    const int i_ub = min(mb_shift + b1, N);

                    const int id_lb = (i_lb == 0) ? 0 : (i_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;
                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int loc_shift = id * db1;
                        const int bd_shift = loc_shift + mb_shift;
                        const int lb = max(lambda, bd_shift) - bd_shift;
                        const int ub = min(db1, i_ub - bd_shift);
                        // Width of block is 'b2',
                        // because memory was allocated for complete matrices
                        double* Ydblock = Y + bd_shift*b2 + lb*dblock_width_WY;
                        double* Wdblock = W + (Ydblock - Y);
                        const double* v_shifted = v + bd_shift + lb;

                        for (int i = 0; i < ub - lb; ++i)
                        {
                            double v_buf = v_shifted[i];
                            *Ydblock = v_buf;
                            *Wdblock = beta * v_buf;
                            Ydblock += dblock_width_WY;
                            Wdblock += dblock_width_WY;
                        }
                    }
                }
            }
            else
            {
                // Computation of production w = (Y^t) * v .
                // Pattern of data access is as follows:
                // 1. Fix big block (ib)
                // 2. Fix small block stripe in fixed big block (id)
                // 3. Fix small block in fixed stripe (jd)
                // 4. Fix row in fixed small block (i) and element of vector 'v'
                // 5. Produce passing of fixed row of Y (j) and vector 'w',
                //    using fixed element of vector 'v'
                for (int ib = b_begin; ib < col_b_count; ++ib)
                {
                    const int mb_shift = ib*b1 + d1_shift;
                    const int ib_lb = max(full_shift, mb_shift) - mb_shift;
                    const int ib_ub = min(mb_shift + b1, N);
                    const int ib_block_h = min(b1, N - mb_shift);
                    const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;

                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int loc_shift = id * db1;
                        const int bd_shift = loc_shift + mb_shift;
                        const int i_lb = max(full_shift, bd_shift) - bd_shift;
                        const int i_ub = min(db1, ib_ub - bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int id_block_h = min(db1, ib_block_h - loc_shift);
                        const double* Ydstripe = Y + bd_shift*b2;
                        const double* v_shifted = v + bd_shift + i_lb;

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it - loc_drow_shift, db2);
                            const int dblock_w = min(db2, b2 - loc_drow_shift);
                            const double* Yrow = Ydstripe +
                                id_block_h * loc_drow_shift + i_lb * dblock_w;
                            double* w_shifted = w + loc_drow_shift;

                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                const double vi = v_shifted[i];
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    w_shifted[j] += Yrow[j] * vi;
                                }
                                Yrow += dblock_w;
                            }
                        }
                    }
                }

                // Computation of production W * ((Y^t) * v) <=> W * w .
                // Pattern of data access is as follows:
                // 1. Fix big block (ib)
                // 2. Fix small block stripe in fixed big block (id)
                // 3. Fix small block in fixed stripe (jd)
                // 4. Fix row in fixed small block (i)
                // 5. Produce passing of fixed row of W (j) and vector 'w'
                const int set_size_SUM_VEC = db1 * sizeof(double);
                for (int ib = 0; ib < col_b_count; ++ib)
                {
                    const int mb_shift = ib*b1 + d1_shift;
                    const int ib_lb = max(lambda, mb_shift) - mb_shift;
                    const int ib_ub = min(mb_shift + b1, N);
                    const int ib_block_h = min(b1, N - mb_shift);

                    const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
                    const int id_ub = (ib == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;
                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int loc_shift = id * db1;
                        const int bd_shift = loc_shift + mb_shift;
                        const int i_lb = max(lambda, bd_shift) - bd_shift;
                        const int i_ub = min(db1, ib_ub - bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int id_block_h = min(db1, ib_block_h - loc_shift);
                        double* Wdstripe = W + bd_shift * b2;
                        double* Ydstripe = Y + (Wdstripe - W);
                        // Vector to contain results of production of W's row and vector w
                        memset(sum_vec, 0, set_size_SUM_VEC);

                        for (int jd = 0; jd < dblock_col_index + 1; ++jd)
                        {
                            const int loc_drow_shift = jd * db2;
                            const int j_ub = min(it, loc_drow_shift + db2) - loc_drow_shift;
                            const int dblock_w = min(db2, b2 - loc_drow_shift);
                            double* Wrow = Wdstripe +
                                id_block_h * loc_drow_shift + i_lb * dblock_w;
                            const double* w_shifted = w + loc_drow_shift;
                            double* sum_vec_shifted = sum_vec + i_lb;

                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                double sum_loc = 0.0;
                                for (int j = 0; j < j_ub; ++j)
                                {
                                    sum_loc += Wrow[j] * w_shifted[j];
                                }
                                // Scalar production of W matrix line and vector w
                                sum_vec_shifted[i] += sum_loc;

                                Wrow += dblock_w;
                            }
                        }
                        // Insertion of computed values into 'it'-th column of Y and W
                        const double* v_shifted = v + bd_shift + i_lb;
                        // Pointers to small blocks containing 'it'-th column of Y and W
                        double* Wdblock_add = Wdstripe + dblock_col_index*id_block_h*db2 +
                            i_lb*dblock_width_WY + loc_it;
                        double* Ydblock_add = Ydstripe + (Wdblock_add - Wdstripe);
                        const double* sum_vec_shifted = sum_vec + i_lb;

                        for (int i = 0; i < i_iter_count; ++i)
                        {
                            double buf = v_shifted[i];
                            *Wdblock_add = beta * (buf + sum_vec_shifted[i]);
                            *Ydblock_add = buf;
                            Wdblock_add += dblock_width_WY;
                            Ydblock_add += dblock_width_WY;
                        }
                    }
                }
            }
        }

        // * Transformation of remaining block-columns of A...
        // A(lamda:M-1, t:N-1) = (I + W*(Y^t))^t * A(lamda:M-1, t:N-1) (==)
        // (==) A(lamda:M-1, t:N-1) + [Y*(W^t)] * A(lamda:M-1, t:N-1)

        set_minor_double_block(WY, 0, N, b1, b1, db1, db1, d1_shift, d1_shift);

        // Multiplication Y * W^t .
        // Remark: WY matrix get (B1,B1,D1,D1)-layout.
#pragma omp parallel for schedule(dynamic, 1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib * b1;
            const int ib_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int ib_ub = min(i_mb_shift + b1, N);
            const int ib_stripe_h = min(b1, N - i_mb_shift);
            double* WYstripe = WY + i_mb_shift*N;

            const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
            const int id_ub = (ib == col_b_count - 1) ?
                dblock_count_in_diff_bcol :
                dblock_count_in_bcol;
            for (int jb = 0; jb < col_b_count; ++jb)
            {
                const int j_mb_shift = d1_shift + jb *b1;
                const int jb_lb = max(lambda, j_mb_shift) - j_mb_shift;
                const int jb_ub = min(j_mb_shift + b1, N);
                const int jb_stripe_h = min(b1, N - j_mb_shift);
                double* WYblock = WYstripe + ib_stripe_h*j_mb_shift;

                const int jd_lb = (jb_lb == 0) ? 0 : (jb_lb / db1);
                const int jd_ub = (jb == col_b_count - 1) ?
                    dblock_count_in_diff_bcol :
                    dblock_count_in_bcol;
                for (int id = id_lb; id < id_ub; ++id)
                {
                    const int i_loc_shift = id * db1;
                    const int i_bd_shift = i_mb_shift + i_loc_shift;
                    const int i_lb = max(lambda, i_bd_shift) - i_bd_shift;
                    const int i_ub = min(db1, ib_ub - i_bd_shift);
                    const int i_iter_count = i_ub - i_lb;
                    const int id_block_h = min(db1, ib_stripe_h - i_loc_shift);
                    const double* Ydstripe = Y + i_bd_shift*b2;
                    double* WYdstripe = WYblock + i_loc_shift*jb_stripe_h;
                    
                    // Upper bound for 'kd'-level small-block-cycle is defined by
                    // lower triangular form of matrix Y.
                    // If 'id'-th small-block-stripe is crossing matrix diagonal,
                    // then upper bound for small-block-stripe passing equals
                    // index of small block, that is containing diagonal element
                    // with maximum index within 'id'-th small-block-stripe height limits.
                    const int kd_ud = (i_bd_shift < t) ?
                        ((min(i_bd_shift + db1, ib_ub) - lambda) / db2 + 1):
                        dblock_count_in_curr_block_row;

                    for (int jd = jd_lb; jd < jd_ub; ++jd)
                    {
                        const int j_loc_shift = jd * db1;
                        const int j_bd_shift = j_mb_shift + j_loc_shift;
                        const int j_lb = max(lambda, j_bd_shift) - j_bd_shift;
                        const int j_ub = min(db1, jb_ub - j_bd_shift);
                        const int j_iter_count = j_ub - j_lb;
                        const int jd_block_h = min(db1, jb_stripe_h - j_loc_shift);
                        const int WYrow_diff_val = jd_block_h - j_iter_count;
                        const int WY_diff_val = i_iter_count * jd_block_h;
                        const double* Wdstripe = W + j_bd_shift*b2;
                        double* WYrow_shifted = WYdstripe + id_block_h*j_loc_shift +
                            i_lb * jd_block_h + j_lb;

                        for (int kd = 0; kd < kd_ud; ++kd)
                        {
                            const int k_loc_shift = kd * db2;
                            const int k_ub = min(d2 - k_loc_shift, db2);
                            const int kd_dblock_w = min(db2, b2 - k_loc_shift);
                            const int Wrow_diff_val = kd_dblock_w - k_ub;
                            const int W_diff_val = j_iter_count * kd_dblock_w;
                            const double* Yrow = Ydstripe +
                                id_block_h * k_loc_shift + i_lb * kd_dblock_w;
                            const double* Wrow = Wdstripe +
                                jd_block_h * k_loc_shift + j_lb*kd_dblock_w;

                            for (int i = 0; i < i_iter_count; ++i)
                            {
                                for (int j = 0; j < j_iter_count; ++j)
                                {
                                    double sum = 0.0;
                                    for (int k = 0; k < k_ub; ++k)
                                    {
                                        sum += (*Yrow++) * (*Wrow++);
                                    }
                                    (*WYrow_shifted++) += sum;

                                    Yrow -= k_ub;
                                    Wrow += Wrow_diff_val;
                                }

                                WYrow_shifted += WYrow_diff_val;
                                Yrow += kd_dblock_w;
                                Wrow -= W_diff_val;
                            }

                            WYrow_shifted -= WY_diff_val;
                        }
                    }
                }
            }
        }

        copy_minor_double_block(Abuf, A, N, b1, b2, db1, db2, d1_shift, d2_shift);

        // Multiplication A(lamda:N-1, t:N-1) += (Y * W^t) * A(lamda:M-1, t:N-1) <=>
        //                A(lamda:N-1, t:N-1) += WY * A(lamda:N-1, t:N-1) .
        // Layout pattern: (b1 x b2, d1 x d2) = (b1 x b1, d1 x d1) * (b1 x b2, d1 x d2).
#pragma omp parallel for schedule(dynamic, 1)
        for (int ib = 0; ib < col_b_count; ++ib)
        {
            const int i_mb_shift = d1_shift + ib*b1;
            const int ib_lb = max(lambda, i_mb_shift) - i_mb_shift;
            const int ib_ub = min(i_mb_shift + b1, N);
            const int ib_stripe_h = min(b1, N - i_mb_shift);
            const double* WY_stripe = WY + i_mb_shift * N;
            double* Abuf_stripe = Abuf + (WY_stripe - WY);

            const int id_lb = (ib_lb == 0) ? 0 : (ib_lb / db1);
            const int id_ub = (ib == col_b_count - 1) ?
                dblock_count_in_diff_bcol :
                dblock_count_in_bcol;
            for (int jb = 0; jb < row_b_count; ++jb)
            {
                const int j_mb_shift = t + jb * b2;
                const int jb_ub = min(j_mb_shift + b2, N);
                const int jb_column_w = min(b2, N - j_mb_shift);
                double* Abuf_block = Abuf_stripe + j_mb_shift * ib_stripe_h;

                const int jd_ub = (jb == row_b_count - 1) ?
                    dblock_count_in_diff_brow :
                    dblock_count_in_brow;
                for (int kb = 0; kb < col_b_count; ++kb)
                {
                    const int k_mb_shift = d1_shift + kb * b1;
                    const int kb_lb = max(lambda, k_mb_shift) - k_mb_shift;
                    const int kb_ub = min(k_mb_shift + b1, N);
                    const int kb_stripe_h = min(b1, N - k_mb_shift);
                    const double* WY_block = WY_stripe + ib_stripe_h * k_mb_shift;
                    const double* A_block = A + k_mb_shift * N +
                        j_mb_shift * kb_stripe_h;

                    const int kd_lb = (kb_lb == 0) ? 0 : (kb_lb / db1);
                    const int kd_ub = (kb == col_b_count - 1) ?
                        dblock_count_in_diff_bcol :
                        dblock_count_in_bcol;
                    for (int id = id_lb; id < id_ub; ++id)
                    {
                        const int i_loc_shift = id * db1;
                        const int i_bd_shift = i_mb_shift + i_loc_shift;
                        const int i_lb = max(lambda, i_bd_shift) - i_bd_shift;
                        const int i_ub = min(db1, ib_ub - i_bd_shift);
                        const int i_iter_count = i_ub - i_lb;
                        const int id_block_h = min(db1, ib_stripe_h - i_loc_shift);
                        const double* WY_dstripe = WY_block + i_loc_shift*kb_stripe_h;
                        double* Abuf_dstripe = Abuf_block + i_loc_shift*jb_column_w;

                        for (int jd = 0; jd < jd_ub; ++jd)
                        {
                            const int j_loc_shift = jd * db2;
                            const int j_bd_shift = j_mb_shift + j_loc_shift;
                            const int j_ub = min(j_bd_shift + db2, jb_ub) - j_bd_shift;
                            const int jd_column_w = min(db2, jb_column_w - j_loc_shift);
                            const int Abuf_diff_val = i_iter_count * jd_column_w;
                            double* Abuf_row = Abuf_dstripe +
                                id_block_h * j_loc_shift + i_lb * jd_column_w;

                            for (int kd = kd_lb; kd < kd_ub; ++kd)
                            {
                                const int k_loc_shift = kd * db1;
                                const int k_bd_shift = k_mb_shift + k_loc_shift;
                                const int k_lb = max(lambda, k_bd_shift) - k_bd_shift;
                                const int k_ub = min(k_bd_shift + db1, kb_ub) - k_bd_shift;
                                const int k_iter_count = k_ub - k_lb;
                                const int A_row_diff_val = k_iter_count * jd_column_w - 1;
                                const int kd_block_h = min(db1, kb_stripe_h - k_loc_shift);
                                const double* WY_row_shifted = WY_dstripe +
                                    id_block_h*k_loc_shift + i_lb*kd_block_h + k_lb;
                                const double* A_row = A_block +
                                    k_loc_shift * jb_column_w +
                                    kd_block_h * j_loc_shift + k_lb * jd_column_w;

                                for (int i = 0; i < i_iter_count; ++i)
                                {
                                    for (int j = 0; j < j_ub; ++j)
                                    {
                                        double sum = 0.0;
                                        for (int k = 0; k < k_iter_count; ++k)
                                        {
                                            sum += WY_row_shifted[k] * (*A_row);
                                            A_row += jd_column_w;
                                        }
                                        Abuf_row[j] += sum;

                                        A_row -= A_row_diff_val;
                                    }

                                    WY_row_shifted += kd_block_h;
                                    Abuf_row += jd_column_w;
                                    A_row -= j_ub;
                                }
                                Abuf_row -= Abuf_diff_val;
                            }
                        }
                    }
                }
            }
        }

        copy_minor_double_block(A, Abuf, N, b1, b2, db1, db2, d1_shift, d2_shift);
    }// WHILE
    
    free(work_memory);

    return A;
}

#ifdef __cplusplus
    }
#endif
