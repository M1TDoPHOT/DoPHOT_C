#ifndef CHISQ_H_
#define CHISQ_H_

/* dophot function converted to c double function 02-21-2012 */

double chisq_( double (*FUNCTN)(short int*, float*, float*, int*, int*), short int** X, float* Y, float* YE, int* N_ptr, float* A, float* FA, float* C_ptr, int* M_ptr, float* ACC, float* ALIM, int* IT_ptr);
/* dimensions
      DIMENSION A(m+1), FA(m+1), C(M+1,M+1), ACC(m+1), ALIM(m+1)
      DIMENSION X(n), Y(n), YE(n)
*/

#endif

