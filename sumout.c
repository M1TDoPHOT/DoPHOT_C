#include <stdio.h>
#include "tuneable.h"
#include "estuff_struct.h"
#include "eoff_struct.h"
#include "string.h"
#include "cast_arr.h"
#include "outputs.h"

/* dophot subroutine converted to c void function 02-06-2012 */

char* sumout_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPARR_ptr, float* APPARR, int* lNAPPLE_ptr, float* PROBG_ptr, int which_model)
{

     /* dereference the pointers */
     int   I      = *I_ptr;
     int   ITYPE  = *ITYPE_ptr;
//unused     int   NPARR  = *NPARR_ptr;
//unused     int   lNAPPLE = *lNAPPLE_ptr;
     float PROBG  = *PROBG_ptr; 
     
     /* renaming used common block variables */
     float XFRAC = eoff_.xfrac;
     float YFRAC = eoff_.yfrac;
     short int* EMSUB = estuff_.emsub;
     float*     EMERR = estuff_.emerr;
     float**    EMPAR = estuff_.empar;

     /* substace of subroutine begins here */
     static char OUTSTRING[300];
     float t1, t2, t3, t4;
     int K = I-1;
     if (EMSUB[K] != 0){
          t1 = EMPAR[K][0];
          t2 = EMPAR[K][1]*10000.0f;
          t3 = EMPAR[K][2] + XFRAC;
          t4 = EMPAR[K][3] + YFRAC;
     }
     else{
          t1 = EMPAR[K][0];
          t2 = EMPAR[K][1];
          t3 = EMPAR[K][2];
          t4 = EMPAR[K][3];
     }

     char** flags = tune16_.flags;

     if ((strncmp(flags[0], "PGAUSS", 5) == 0) ||
         (strncmp(flags[0], "GAUSS",  5) == 0)) {
          sprintf(OUTSTRING," %4d %2d %10.4E %10.2f %10.2f %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %7.3f %7.3f %10.2E %10.3E %10.3E %10.2f %10.2f %7.3f\n",
          I, ITYPE,  STPARR[0], STPARR[2], STPARR[3],
          STPARR[1], STPARR[4], STPARR[5], STPARR[6], 
          APPARR[0], APPARR[1], APPARR[2], APPARR[3],
          PROBG, t1, t2, t3, t4, EMERR[K]);
     }
     if ((strncmp(flags[0], "PGALTMODEL", 5) == 0) ||
         (strncmp(flags[0], "GALTMODEL",  5) == 0)) {
          sprintf(OUTSTRING," %4d %2d %10.4E %10.2f %10.2f %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %2d %10.3E %10.3E %7.3f %7.3f %10.2E %10.3E %10.3E %10.2f %10.2f %7.3f\n",
          I, ITYPE,  STPARR[0], STPARR[2], STPARR[3],
          STPARR[1], STPARR[4], STPARR[5], STPARR[6], 
          STPARR[7], STPARR[8], STPARR[9], STPARR[10],
          which_model,
          APPARR[0], APPARR[1], APPARR[2], APPARR[3],
          PROBG, t1, t2, t3, t4, EMERR[K]);
     }
     if( strncmp(flags[0], "SERSIC", 5) == 0)  {
          sprintf(OUTSTRING," %4d %2d %10.4E %10.2f %10.2f %10.3E %10.3E %10.3E %10.3E %10.2f %2d %10.3E %10.3E %7.3f %7.3f %10.2E %10.3E %10.3E %10.2f %10.2f %7.3f\n",
          I, ITYPE,  STPARR[0], STPARR[2], STPARR[3],
          STPARR[1], STPARR[4], STPARR[5], STPARR[6], 
          STPARR[7],
          which_model,
          APPARR[0], APPARR[1], APPARR[2], APPARR[3],
          PROBG, t1, t2, t3, t4, EMERR[K]);
     }
     if( strncmp(flags[0], "EXTPGAUSS", 5) == 0)  {
          sprintf(OUTSTRING," %4d %2d %10.4E %10.2f %10.2f %10.3E %10.3E %10.3E %10.3E %10.2f %10.2f %2d %10.3E %10.3E %7.3f %7.3f %10.2E %10.3E %10.3E %10.2f %10.2f %7.3f\n",
          I, ITYPE,  STPARR[0], STPARR[2], STPARR[3],
          STPARR[1], STPARR[4], STPARR[5], STPARR[6], 
          STPARR[7], STPARR[8], 
          which_model,
          APPARR[0], APPARR[1], APPARR[2], APPARR[3],
          PROBG, t1, t2, t3, t4, EMERR[K]);
     }


     return OUTSTRING;
}


