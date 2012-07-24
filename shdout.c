#include <stdio.h>
#include "tuneable.h"
#include "string.h"
#include "outputs.h"

/* dophot subroutine converted to c void function 02-06-2012 */

char* shdout_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPARR_ptr, int which_model)
{

     /* dereference the pointers */
     int   I      = *I_ptr;
     int   ITYPE  = *ITYPE_ptr;
//unused     int   NPARR  = *NPARR_ptr;
     
     /* substace of subroutine begins here */
     static char OUTSTRING[300];
     char** flags = tune16_.flags;

     if ((strncmp(flags[0], "PGAUSS", 5) == 0) ||
         (strncmp(flags[0], "GAUSS",  5) == 0)) {
          sprintf(OUTSTRING," %4d %2d %10.4E %10.2f %10.2f %10.3E %10.3E %10.3E %10.3E\n",
          I, ITYPE,  STPARR[0], STPARR[2], STPARR[3],
          STPARR[1], STPARR[4], STPARR[5], STPARR[6]); 
     }
     if ((strncmp(flags[0], "PGALTMODEL", 5) == 0) ||
         (strncmp(flags[0], "GALTMODEL",  5) == 0)) {
          sprintf(OUTSTRING," %4d %2d %10.4E %10.2f %10.2f %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %2d\n",
          I, ITYPE,  STPARR[0], STPARR[2], STPARR[3],
          STPARR[1], STPARR[4], STPARR[5], STPARR[6], 
          STPARR[7], STPARR[8], STPARR[9], STPARR[10], which_model); 
     }
     if (strncmp(flags[0], "SERSIC", 5) == 0){ 
          sprintf(OUTSTRING," %4d %2d %10.4E %10.2f %10.2f %10.3E %10.3E %10.3E %10.3E %10.2f %2d\n",
          I, ITYPE,  STPARR[0], STPARR[2], STPARR[3],
          STPARR[1], STPARR[4], STPARR[5], STPARR[6], 
          STPARR[7], which_model); 
     }
     if (strncmp(flags[0], "EXTPGAUSS", 5) == 0){ 
          sprintf(OUTSTRING," %4d %2d %10.4E %10.2f %10.2f %10.3E %10.3E %10.3E %10.3E %10.2f %10.2f %2d\n",
          I, ITYPE,  STPARR[0], STPARR[2], STPARR[3],
          STPARR[1], STPARR[4], STPARR[5], STPARR[6], 
          STPARR[7], STPARR[8], which_model); 
     }

     return OUTSTRING;
}


