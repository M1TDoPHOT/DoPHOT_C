#ifndef DIAG_MMULT_H_
#define DIAG_MMULT_H_

//routines to multiply one matrix by a diagonal matrix
// on either the left or right or both sides

// multiply matrix mat by a diagonal matrix on the left
void left_diag_mmult(float** mat, float* d_mat, float** result, int dimen);

// multiply matrix mat by a diagonal matrix on the right 
void right_diag_mmult(float** mat, float* d_mat, float** result, int dimen);

// multiply matrix mat by a diagonal matrix on both sides
void rightleft_diag_mmult(float** mat, float* d_mat, float** result, int dimen);

#endif
