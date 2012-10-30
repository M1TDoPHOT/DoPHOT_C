#include <stdlib.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "fitting_matrices_struct.h"
#include "mini_mathlib.h"
#include "cast_arr.h"
#include "guess.h"
#include "lu_comp.h"
#include "elarea.h"
#include "probgal.h"

/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:I THINK I WANT TO CALCULATE THE DERIVATIVE OF THE BEST FITTING CHISQUARED
C:WITH RESPECT TO A VARIATION IN A(2),A(5),A(6) AND A(7) WHICH KEEPS FLUX
C:CONSTANT.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

/* dophot function converted to c double function 02-20-2012 */
//onestar is a pointer to the function name ONESTAR which returns double
//onestar is usually pseud2d

double probgal_(double (*ONESTAR)(short int*, float*, float*, int*, int*), short int** XX, float* YY, float* ZZ, int* NPT_ptr, float* A, float* FA)
{
     /* dereferencing the pointer */
     int NPT = *NPT_ptr;
     /* renaming used common block variable */
     int NPAR  = tune4_.npar;
     int lverb = tune14_.lverb;

     /* substance of fucntion begins here */
     float** C  = fitting_matrices_.c_mat;
     float** LU = fitting_matrices_.lu;
     float* V   = fitting_matrices_.v;
     float* Vsol = fitting_matrices_.vsol;
     int* INDX = fitting_matrices_.index;
     float* C2 = fitting_matrices_.c_list; //NPMAX
     // none gaurenteed to be 0'd

     int I,J,K; 
     float CHI1 = 0.0f;
     double F_double;
     float F, FAJ;
     int NPAR_dum  = NPAR; //needed for passes to ONESTAR
     int ZERO_dum  = 0; //needed for passes to ONESTAR
     float DCHI = 0.0f;
     float B4, B5, B6;
     float AAREA, BAREA;
     float CHIPER; 
     float PROBGAL; 
    
     //zero C and C2
     for(I = 0; I < NPMAX; I++){
          C2[I] = 0.0f;
          for(J = 0; J < NPMAX; J++){
               C[I][J] = 0.0f;
          } 
     } 

     for(I = 0; I < NPT; I++){
          F_double = (*ONESTAR)(XX[I], A, FA, &NPAR_dum, &ZERO_dum);
          F = (float)(F_double) - YY[I];
          CHI1 += F*F/ZZ[I];
          for(J = 0; J < NPAR; J++){
               FAJ = FA[J]/ZZ[I];
               C2[J] += FAJ*F;
               V[J] = C2[J];
               for(K = 0; K <= J; K++){
                    C[K][J] += FAJ*FA[K];
               }
          }//end J loop
     }//end I loop

     /* diagonalizing C*/
     for (I = 0; I < NPAR; I++){
          for (J = 0; J <= I; J++){
               C[I][J] = C[J][I];
          }
     }
     /* matrix LU decomposition and solving */
     lu_decompose(C, NPAR, INDX, LU); 
     lu_solve(LU, NPAR, INDX, V, Vsol); 
     for (J = 0; J < NPAR; J++){
          V[J]    = Vsol[J];
          Vsol[J] = 0.0f;
     }

     for (I = 0; I < NPAR; I++){
          DCHI -= C2[I]*V[I];
          for (J = 0; J < NPAR; J++){
               DCHI += V[I]*C[I][J]*V[J]/2.0f;
          }
     }
     B4 = A[4] - V[4];
     B5 = A[5] - V[5];
     B6 = A[6] - V[6];
     AAREA = elarea_(A+4, A+5, A+6);
     BAREA = elarea_(&B4, &B5, &B6);
     if (DCHI > 0.0f){
          if (lverb > 20){
               fprintf(logfile, "Object at %f %f has -ve Delta Chi**2 !\n",
                                         A[2], A[3]);
               fprintf(logfile, "DCHI = %f \n", DCHI);
          }
          PROBGAL = 0.0f;
     }
     else{
          if (AAREA <= 0.0f){
               PROBGAL = 0.0f;
          }
          else{
               CHIPER = CHI1/NPT;
               PROBGAL = fabsf(DCHI/CHIPER) * (BAREA-AAREA)/fabsf(BAREA-AAREA);
          }
     }

     return PROBGAL;
}
          



	
