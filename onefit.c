#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "cast_arr.h"
#include "chisq.h"
#include "diag_mmult.h"
#include "onefit.h"

/* a wrapper function for chisq that will convert parameters and limits to
 and from log space if desired for fitting purposes */
// currently, update for each model switch, latter make toggle with flag[0]

double onefit_( double (*FUNCTN)(short int*, float*, float*, int*, int*), short int** X, float* Y, float* YE, int* N_ptr, float* A, float* FA, float* C_ptr, int* M_ptr, float* ACC, float* ALIM, int* IT_ptr)
{
     int M = *M_ptr; //number of fit parameters
     float** C    = malloc_float_2darr(M,M);
     float* D_ptr = malloc_float_1darr(M*M);
     float** D    = malloc_float_2darr(M,M);

     float* J     = malloc_float_1darr(M); //Jacobian, assuming diagonal
                        // and only keeping diagonal elements

     int i;
     float ALIM1, ACC1;
     double chisq_return;

     /* move intensity to log space for fitting
        linear space variable I.  log space variable N.
        currently, A[1]  = Intensity I = exp(N);
        currently, FA[1] = 1/denom */

     // set value
     if( A[1] > 0.00f ){
          A[1] = logf(A[1]);
     }
     else{
          fprintf(logfile, 
               "caution: negative intensity in fit, %f\n", A[1]);
          A[1] = -9.999f;
     } // now A[1] = N = log(I).  I = exp(N) 

     // set Derivative
     FA[1] = A[1]*FA[1]; //now FA[1] = Intensity/denom = exp(N)/denom

     // set limits
     ACC1  = ACC[1];
     ALIM1 = ALIM[1];
     ACC[1]  = -0.01f;
     ALIM[1] = -10.0f;

     // C_ptr is updated by chisq_ not used by it, 
     // so no need to set covariance, only reassign it
     
     /* end move intensity to log space for fitting */

     /* do fit */
     chisq_return = chisq_(FUNCTN, X, Y, YE, N_ptr, A, FA, D_ptr, 
                          M_ptr, ACC, ALIM, IT_ptr);

     /* move intensity back to linear space
        currently A[1] = N = log(I).  I = exp(N);
        currently FA[1] = I/denom */

     // set value
     A[1]  = expf(A[1]); // now A  = Intensity I 

     // set Derivative
     FA[1] = FA[1]/A[1]; // now FA = 1/denom

     // reset limits
     ACC[1]  = ACC1;
     ALIM[1] = ALIM1;

     // set Jacobian
     for(i = 0; i < M; i++){
          J[i] = 1.0f;
     }
     J[1] = A[1]; 

     // set covariance
     recast_float_1dto2darr(M, M, D_ptr, D);
     rightleft_diag_mmult(D, J, C, M);
     C[1][0] = C[1][0] / A[1]; // I have no idea why these are necessary
     C[0][1] = C[0][1] / A[1]; // but they are
     C[1][1] = C[1][1] / A[1]; // in order to get correct covariance matrix
     recast_float_2dto1darr(M, M, C_ptr, C);
     
     /* end move intensity back to linear space */

     free(J);
     free(D_ptr);
     free_float_2darr(M, D);
     free_float_2darr(M, C);

     return chisq_return;
}
