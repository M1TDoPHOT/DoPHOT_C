#include <stdlib.h>
#include "logh.h"
#include "tuneable.h"
#include "starlist_struct.h"
#include "search_struct.h"
#include "model_struct.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "paravg.h"

/* dophot subroutine converted to c void function 02-02-2012 */
/* RUN THROUGH STAR LIST AND COMPUTE AVERAGE VALUES FOR SKY, ETC. */

void paravg_()
{
 
     /* rename and cast used common block variables */
     int  NPAR   = tune4_.npar;
     int  lverb  = tune14_.lverb;
     float* AVA  = tune15_.ava; //FIX ME changed in this routine!!!
     int  NSTOT  = search_.nstot;
     int* IMTYPE = starlist_.IMTYPE;
     float** SHADOW  = starlist_.shadow;
     float** SHADERR = starlist_.shaderr;
     int* WHICH_MODEL = model_.which_model;
 
     /* substance of subroutine begins here */ 
     float** SUM = malloc_float_2darr(NPAR, 2);
     int I,J;
     int N = 0;
     int CONDIT;

     for (J = 0; J < NPAR; J++){
          SUM[J][0] = 0.0f;
          SUM[J][1] = 0.0f;
     }

     if (NSTOT >= 1){
          for (I = 0; I < NSTOT; I++){
               CONDIT = ((IMTYPE[I] == 1) || (IMTYPE[I] == 11));
               CONDIT = ((CONDIT) && (WHICH_MODEL[I] == 0));
               if (CONDIT){
                    N += 1;
                    for (J = 0; J < NPAR; J++){
                         if (SHADERR[I][J] == 0.0f){
                              fprintf(logfile,"shaderr = 0, i = %7d\n", I+1);
                              fprintf(logfile,"\n");
                         }
                         SUM[J][0] += 1.0f/SHADERR[I][J];
                         SUM[J][1] += SHADOW[I][J]/SHADERR[I][J];
                    } // end of J loop
               }
          } // end of I loop
          if (N > 0){
               for (J = 0; J < NPAR; J++){
                    AVA[J] = SUM[J][1]/SUM[J][0];
               } 
          }
          if (lverb > 10){
               fprintf(logfile, "AVERAGE values so far of NPAR parameters for stars:\n");   
               for (J = 0; J < NPAR; J++){
                    fprintf(logfile, "%9.5f  ", AVA[J]);
               }
               fprintf(logfile, "\n");   
          }
     }

     /* recasting changed pointers and freeing allocated memory */	
     free_float_2darr(NPAR, SUM);
    
} 
