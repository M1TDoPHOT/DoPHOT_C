#include <stdlib.h>
#include "cast_arr.h"
#include "diag_mmult.h"

//routines to multiply one matrix by a diagonal matrix
// on either the left or right or both sides

// multiply matrix mat by a diagonal 'matrix' on the left
// only pass diagonal elements of diag matrix
void left_diag_mmult(float** mat, float* d_mat, float** result, int dimen)
{
     int n, m;
     for( m = 0; m < dimen; m++){
          for( n = 0; n < dimen; n++){
               result[m][n] = d_mat[m]*mat[m][n];
          }
     }
}

// multiply matrix mat by a diagonal matrix on the right 
// only pass diagonal elements of diag matrix
void right_diag_mmult(float** mat, float* d_mat, float** result, int dimen)
{
     int n, m;
     for( m = 0; m < dimen; m++){
          for( n = 0; n < dimen; n++){
               result[m][n] = d_mat[n]*mat[m][n];
          }
     }
}

// multiply matrix mat by a diagonal matrix on both sides
// only pass diagonal elements of diag matrix
void rightleft_diag_mmult(float** mat, float* d_mat, float** result, int dimen)
{
     int n, m;
     for( m = 0; m < dimen; m++){
          for( n = 0; n < dimen; n++){
               result[m][n] = d_mat[n]*d_mat[m]*mat[m][n];
          }
     }
}
