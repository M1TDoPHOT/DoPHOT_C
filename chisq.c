#include <stdlib.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "fitting_matrices_struct.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "lu_comp.h"
#include "chisq.h"

/* dophot function converted to c double function 02-21-2012 */
// imitation of Bevington's (1969) implementation of Marquardt's algorithm.
// similar to a Levenberg-Marquardt algorithm

double chisq_( double (*FUNCTN)(short int*, float*, float*, int*, int*), short int** X, float* Y, float* YE, int* N_ptr, float* A, float* FA, float* C_ptr, int* M_ptr, float* ACC, float* ALIM, int* IT_ptr)
/* dimensions
      DIMENSION A(M), FA(M), C(M,M), ACC(M), ALIM(M)  M=number params
      DIMENSION X(N, 2), Y(N), YE(N) N=number stars
*/
{
     /* dereferencing pointers and 2x2 array pointers */
     int N  = *N_ptr ;
     int M  = *M_ptr ;
     int IT = *IT_ptr;
     float** C  = fitting_matrices_.c_mat;
     recast_float_1dto2darr(M, M, C_ptr, C);
     
     /* renaming used common block variables */
     int lverb   = tune14_.lverb;

     /* substance of subroutine beigns here */
     if (M > MMAX-1){
          printf("chisq thinks M>MMAX-1... possible bug\n");
          return 0;
     }
     
     int CONV, MARQ, LIMIT, skip_ahead;
     float** LU   = fitting_matrices_.lu;
     float* V     = fitting_matrices_.v;
     float* Vsol  = fitting_matrices_.vsol;
     float** B = fitting_matrices_.b_mat;
     int* INDX = fitting_matrices_.index;
     int ZERO_dum = 0; //needed only for function pass
     int FITCALL = 1; //calling the fitting function from chisq
     int J, KK, L, count;
     int IFACT, I, K;
     float CHI, CHIOLD;
     double functn_return;
     float F, F2, YE1, FACT, FAKK;
     float PERDEG, TEMP;
     float CHISQ;

     CONV  = 0; //false
     LIMIT = 0; //false

     for( J = 0; J < M; J++){
          /* CHECK IF PARAMETERS EXCEEDS LIMITS AT THE START.  CHECK ONLY
             ABSOLUTE PARAMETERS (ALIM(J) < 0). */
          if (ALIM[J] < 0.0f){
               LIMIT = (fabsf(A[J]) > fabsf(ALIM[J]));
               if ((LIMIT) && (lverb > 20)){
                    fprintf(logfile,"SELF-DECEPTION HAS OCCURED: INITIAL LIMITS.\n");
               }
          }
     }
     if (lverb > 30){
          fprintf(logfile,"%3d %2d ", 0, 0);
          for(KK = 0; KK < M; KK++){
               fprintf(logfile,"%10.6e ", A[KK]);
          }
          fprintf(logfile,"\n");
     }

     IFACT = 0;
     I     = 0;
     skip_ahead = 0; //false
     while( (I < IT) && (!CONV) && (!LIMIT) && (!skip_ahead) ){
          CHI = 0.0f;

          //zero B array except for B[?][M]
          for(J = 0; J < M; J++){
               B[M][J] = 0.0f;
               for(KK = 0; KK < M; KK++){
                    B[KK][J] = 0.0f;
               } 
          }//end J loop

          for(J = 0; J < N; J++){
               //straightforward SUM (model - data)^2/error
               functn_return = (*FUNCTN)(X[J], A, FA, &M, &FITCALL);
               F = (float)functn_return - Y[J];
               F2 = F*F;
               YE1 = 1.0f/YE[J];
               CHI += F2*YE1;

               //Jacobian transpose * Jacobian /err^2 stored in upper half of 
               //B[0-(M-1)][0-(M-1)],
               //equivalently SUM derivative_a*derivative_b/error in space [a][b]
               //B[M][0-(M-1)] stores Jacobian transpose* (function - val)/err, 
               //equivalently SUM derivative*function/err
               // where sums are over all data points (pixels)
               for(KK = 0; KK < M; KK++){
                    if (FA[KK] != 0.0f){
                         FAKK = FA[KK]/YE[J];
                         B[M][KK] += FAKK*F; 
                         for( L = 0; L <= KK; L++){
                              if (FA[L] != 0.0f){
                                   B[L][KK] += FAKK*FA[L];
                              }
                         } //end L loop
                    }
               }// end KK loop

          }//end J loop

          CHIOLD = CHI;
          K = 0;
          MARQ = 0; //false
          while( (K < 10) && (!MARQ) && (!LIMIT) && (!skip_ahead) ){
               CONV = (K == 0);
               if (CONV){
                    FACT = 0.0f;
               }
               else{
                    FACT = powf(2.0, IFACT);
               }

               // set C to B, with scaling if needed
               // make C full diagonal matrix rather than justt upper half of
               // Jacobian transpose Jacobian
               for (J = 0; J < M; J++){
                    for(L = 0; L < J; L++){ //look here
                         C[L][J] = B[L][J];
                         C[J][L] = B[L][J];
                    }
                    C[J][J] = (1+FACT)*B[J][J];
                    V[J]    = B[M][J];
               } // end J loop

               /* matrix LU decomposition and solving 
                  for Cx = V, replace V with solution x when finished*/
               lu_decompose(C, M, INDX, LU);
               lu_solve(LU, M, INDX, V, Vsol);
               for (J = 0; J < M; J++){
                    V[J]    = Vsol[J];
                    Vsol[J] = 0.0f;
               }
               /* V now contains best guess delta for parameters */

               for (J = 0; J < M; J++){
                    A[J] -= V[J];
                    /* CHECK IF CHANGE IN PARAMETERS EXCEEDS LIMITS.  
                       IF ALIM(J) > 0, THEN CONSIDER FRACTIONAL CHANGES.  
                       IF ALIM(J) < 0, CONSIDER ABSOLUTE CHANGES.  
                       IF ALIM(J) = 0, IGNORE THIS TEST. */
                    if (ALIM[J] > 0.0f){
                         LIMIT = (LIMIT || (fabsf(V[J]/A[J]) > ALIM[J]) );
                         if (LIMIT && (lverb > 20)){
                              fprintf(logfile, "SELF-DECEPTION HAS OCCRED: FRAC LIMITS.\n");
                         }
                    }
                    if (ALIM[J] < 0.0f){
                         LIMIT = (LIMIT || (fabsf(A[J]) > fabsf(ALIM[J])) );
                         if (LIMIT && (lverb > 20)){
                              fprintf(logfile, "SELF-DECEPTION HAS OCCRED: ABS LIMITS.\n");
                         }
                    }
                         
                    if (ACC[J] >= 0.0f){
                         CONV = (CONV && (fabsf(V[J]/A[J]) <= ACC[J]) );
                    }
                    else{ 
                         CONV = (CONV && (fabsf(V[J]) <= fabsf(ACC[J])) );
                    }
               } //end J loop

               if (lverb > 30){
                    fprintf(logfile,"%3d %2d ", I+1, K+1);
                    for(KK = 0; KK < M; KK++){
                         fprintf(logfile,"%10.6e ", A[KK]);
                    }
                    fprintf(logfile,"\n");
               }

               if (CONV){
                    MARQ = 1; //true
               }
               else{
                    if (!LIMIT){
                         //compute straightforward chi again, but with new, trial A params
                         CHI = 0.0f;
                         for(J = 0; J < N; J++){
                              functn_return = (*FUNCTN)(X[J], A, FA, &ZERO_dum, &FITCALL);
                              F = (float)functn_return - Y[J];
                              F2 = F*F;
                              YE1 = 1.0f/YE[J];
                              CHI += F2*YE1;
                         }
                         if (K == 1){ 
                              IFACT -= 1;
                         }

                         if (CHI < (1.006f*CHIOLD)){ //keep new params
                              MARQ = 1; //true
                         }
                         else{
                              if (K == 1){ //less trust in derivatives
                                   IFACT += 1;
                              }
                              if (K >= 1){ //even less
                                   IFACT += 1;
                              }
                              if (IFACT > 10){ //no trust
                                   skip_ahead = 1; //true
                              }
                              if (!skip_ahead){ //don't keep current update to A
                                   for(J = 0; J < M; J++){
                                        A[J] += V[J];
                                   }
                              }
                         } //end if CHI if/else
                    } 
               } //end if CONV if/else

               K += 1; //10 attempts to converge
          } // end (K < 10 while)

          I += 1;
     } // end (I < IT while)

     // if the fit is not bad because a limit has been reached in the params
     // compute the true covariance matrix B
     if (!LIMIT){
          //set B to the identity matrix
          for (J = 0; J < M; J++){
               for (I = 0; I < M; I++){
                    B[J][I] = 0.0f;
               }
               B[J][J] = 1.0f;
          }

          for (J = 0; J < M; J++){
               /* matrix solving with LU solution */
               lu_solve(LU, M, INDX, B[J], Vsol);
               for(count = 0; count < M; count++){
                    B[J][count] = Vsol[count];
                    Vsol[count] = 0.0f;
               }
          }
     }


     // make a modified covariance matrix containing the sqrt of the variances 
     // and the normalized correlations.  This is the working 'covariance matrix' C_ptr 
     // use V for storing sqrt(variances)
     if ( CONV && (!LIMIT) ){
          TEMP = (float)(max((N-M), 1));
          PERDEG = sqrtf(CHI/TEMP);

          //copy B to C
          for(I = 0; I < M; I++){
               for(J = 0; J < M; J++){
                    C[I][J] = B[I][J]; 
               }
          }

          //set V to sqrt variances
          for(I = 0; I < M; I++){
               if (B[I][I] > 0.0f){
                    V[I] = fabsf(sqrtf(B[I][I]));
               }
               else{
                    V[I] = 1.0E10f;
                    if (lverb > 20){
                         fprintf(logfile,"TROUBLE: NEGATIVE AUTOVARIANCE FOR I, B(I,I) = %d %f\n",I+1,B[I][I]);
                    }
               }
          }

          // modify the covariance matrix
          for(I = 0; I < M; I++){
               for(J = 0; J < M; J++){
                    C[J][I] = C[J][I]/(V[I]*V[J]); 
               }
               C[I][I] = V[I]*PERDEG;
          } //end I loop
          CHISQ = CHIOLD;
     }
     else{
          CHISQ = 1.0E10f;
     }

     /* recasting changed pointers*/
     recast_float_2dto1darr(M, M, C_ptr, C);

     fprintf(logfile, "chisq returns: %10.6g \n", CHISQ);
     return CHISQ;
}
