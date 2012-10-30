#include <stdlib.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "free_parking_struct.h"
#include "parpred_struct.h"
#include "byvirtue_struct.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "parinterp.h"
#include "galaxy.h"

/* dophot logical function converted to c int function 02-03-2012 */

int galaxy_(float* A, float* ERR, float* STARPAR)
{

     /* renaming used common block variables */
     float* CHI   = byvirtue_.chi;
     float* PARMS = parpred_.parms;
     float* SIG   = tune8_.sig;
     float CHICRIT = tune9_.chicrit;
     int lverb = tune14_.lverb;


     /* substance of function begins here */
     float* B   = free_parking_.npmaxarray_1;
     float* TOT = free_parking_.npmaxarray_2; //though only need 3 elements
     int I;
     float dum, NSIGMA;
     int GALAXY = 0; //FALSE

     for (I = 0; I < 4; I++){
          CHI[I] = 0.0f;
     }
     
     if (   (A[4] >= 2.0f*sqrtf(ERR[4])) 
         && (A[6] >= 2.0f*sqrtf(ERR[6])) 
         && (A[1] >= 2.0f*sqrtf(ERR[1])) ){
          dum = parinterp_(STARPAR+2, STARPAR+3, B);
          TOT[0] = max1(PARMS[3], (SIG[0]*B[4])*(SIG[0]*B[4])); 
          TOT[1] = max1(PARMS[4], (SIG[1]*SIG[1])/(B[4]*B[6]));
          TOT[2] = max1(PARMS[5], (SIG[2]*B[6])*(SIG[2]*B[6])); 
          CHI[0] = (A[4] - B[4])*(A[4] - B[4]) / (TOT[0] + ERR[4]);
          CHI[1] = (A[5] - B[5])*(A[5] - B[5]) / (TOT[1] + ERR[5]);
          CHI[2] = (A[6] - B[6])*(A[6] - B[6]) / (TOT[2] + ERR[6]);
          if (A[4] < B[4]){
               CHI[0] = 0.0f;
          }
          if (A[6] < B[6]){
               CHI[2] = 0.0f;
          }
          CHI[3] = CHI[0] + CHI[1] + CHI[2];
          if (lverb > 20){
               fprintf(logfile,"(GALAXY...) Object at: %11.6f  %11.6f  \n",
                       STARPAR[1], STARPAR[2]);  
               fprintf(logfile,"Chisq-X, Chisq-XY, Chisq-Y, Chisq-TOT: "); 
               for (I = 0; I < 4; I++){
                    fprintf(logfile,"%11.6g  ", CHI[I]);
               }
               fprintf(logfile,"\n");
          }
          NSIGMA = sqrtf(CHI[3]);
          if (lverb > 20){
               fprintf(logfile,"NSIGMA = %9.6g \n", NSIGMA);
          }
          GALAXY = (CHI[3] >= CHICRIT);
     }

     return GALAXY;

}
               


	  
