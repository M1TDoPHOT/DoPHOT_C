#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
// NOTE: the working covariance has the square root of the autovariance along the diagonal, and the correlation on the off-diagonal

double onefit_( double (*FUNCTN)(short int*, float*, float*, int*, int*), short int** X, float* Y, float* YE, int* N_ptr, float* A, float* FA, float* C_ptr, int* M_ptr, float* ACC, float* ALIM, int* IT_ptr, int which_model)
{
     int M = *M_ptr; //number of fit parameters
     float** C    = malloc_float_2darr(M,M);
     float* D_ptr = malloc_float_1darr(M*M);
     float** D    = malloc_float_2darr(M,M);

     float* J     = malloc_float_1darr(M); //Jacobian, assuming diagonal
                        // and only keeping diagonal elements
     float* sqrt_variance = malloc_float_1darr(M); 

     int i, jj;
     float ALIM1, ACC1;
     float ALIM7, ACC7, ALIM8, ACC8;
     double chisq_return;

     int model = 0;  // default gauss, pseudogauss, or not a 
                     // shape fit is 0
     if ((strncmp(tune16_.flags[0], "EXTPGAUSS", 5) == 0) &&
         (which_model == 0)){ //not reverted to pgauss
          model = 1;
     }

     /* move intensity to log space for fitting
        linear space variable I.  log space variable N.
        currently, A[1]  = Intensity I = exp(N);
        currently, FA[1] = 1/denom */

     // set value for A[1]
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

     // set value, derivative and limits for A[7] for extended 
     // pseudogaussian model to log space for fitting
     /* if expanded pseudogaussian model 
        linear space variable beta4.  log space variable R4.
        linear space variable beta6.  log space variable R6.
        currently, A[7]  = beta4 = exp(R4)
        currently, FA[7] = val.  in log space should be beta4*val */ 

     if (model == 1){
          //set value
          if( A[7] > 0.00f ){
               A[7] = logf(A[7]);
          }
          else{
               fprintf(logfile, 
                    "caution: negative beta4 in fit, %f\n", A[7]);
               A[7] = -9.999f;
          } // now A[7] = R4 = log(Beta4).  Beta4 = exp(R4) 
          if( A[8] > 0.00f ){
               A[8] = logf(A[8]);
          }
          else{
               fprintf(logfile, 
                    "caution: negative beta6 in fit, %f\n", A[8]);
               A[8] = -9.999f;
          } // now A[8] = R6 = log(Beta6).  Beta6 = exp(R6) 

          // set Derivative
          FA[7] = A[7]*FA[7]; //now FA[7] = Beta4*val = exp(R4)*val
          FA[8] = A[8]*FA[8]; //now FA[8] = Beta6*val = exp(R6)*val

          // set limits
          ACC7  = ACC[7];
          ALIM7 = ALIM[7];
          ACC[7]  = 0.02f;
          ALIM[7] = -10.0f;
          ACC8  = ACC[8];
          ALIM8 = ALIM[8];
          ACC[8]  = 0.02f;
          ALIM[8] = -10.0f;
     }
          
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
     
     // set value, derivative and limits for A[7] for extended 
     // pseudogaussian model back to linear space
     /* if expanded pseudogaussian model 
        linear space variable beta4.  log space variable R.
        currently, A[7]  = R4 = log(Beta4)
        currently, FA[7] = Beta4*val.  in linear space should be val  
        currently, A[8]  = R6 = log(Beta6)
        currently, FA[8] = Beta6*val.  in linear space should be val */ 
        
     if (model == 1){
          // set value
          A[7]  = expf(A[7]); // now A  = Beta4
          A[8]  = expf(A[8]); // now A  = Beta6
          // if problem in fit, throw back to pgauss values
          if ((isnan(A[7])) || (isinf(A[7]))){
               A[7] = 1.0f;
          }
          if ((isnan(A[8])) || (isinf(A[8]))){
               A[8] = 1.0f;
          }
          // ok because chisq already returned bad value

          // set Derivative
          // if problem in fit, flag derivative as very high
          if ((isnan(FA[7])) || (isinf(FA[7]))){
               FA[7] = 10.0f;
          }
          if ((isnan(FA[8])) || (isinf(FA[8]))){
               FA[8] = 10.0f;
          }

          if (A[7] < 0.0001f){
               FA[7] = FA[7] * 1000000.0;
          } //if prevents NANs in the case of small A[7]
          else{
               FA[7] = FA[7]/A[7]; // now FA = val
          }
          if (A[8] < 0.0001f){
               FA[8] = FA[8] * 1000000.0;
          } //if prevents NANs in the case of small A[8]
          else{
               FA[8] = FA[8]/A[8]; // now FA = val
          }

          // reset limits
          ACC[7]  = ACC7;
          ALIM[7] = ALIM7;
          ACC[8]  = ACC8;
          ALIM[8] = ALIM8;
     }


     // set Jacobian
     for(i = 0; i < M; i++){
          J[i] = 1.0f;
     }
     J[1] = A[1]; 
     if ((model == 1) && (M == tune4_.nfit2)){
          J[7] = A[7]; 
          J[8] = A[8]; 
     }

     // set covariance
     /* NOTE: the working 'covariance' has the square root of the autovariance 
        along the diagonal, and normalized correlations on off axis elements.
        See chisqr.c for details.  */
     recast_float_1dto2darr(M, M, D_ptr, D);
     //set sqrt_variances
     for (i = 0; i < M; i++){
          sqrt_variance[i] = D[i][i];
     }
     //make covariance matrix true covariance matrix
     for (i = 0; i < M; i++){
          D[i][i] = 1.0f;
          for (jj = 0; jj < M; jj++){
               D[i][jj] = D[i][jj]*sqrt_variance[i]*sqrt_variance[jj];
          }
     }
     rightleft_diag_mmult(D, J, C, M);
     //set sqrt_variances
     for (i = 0; i < M; i++){
          sqrt_variance[i] = fabsf(sqrtf(C[i][i]));
     }
     //make true covariance matrix working covariance matrix
     for (i = 0; i < M; i++){
          for (jj = 0; jj < M; jj++){
               C[i][jj] = C[i][jj]/(sqrt_variance[i]*sqrt_variance[jj]);
          }
          C[i][i] = sqrt_variance[i];
     }

     recast_float_2dto1darr(M, M, C_ptr, C);
     
     /* end move intensity back to linear space */

     free(J);
     free(sqrt_variance);
     free(D_ptr);
     free_float_2darr(M, D);
     free_float_2darr(M, C);

     return chisq_return;
}


