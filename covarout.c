#include <stdio.h>
#include "tuneable.h"
#include "string.h"
#include "outputs.h"

char* covarout_(float* matrix)
{

     /* substace of subroutine begins here */
     static char OUTSTRING[300];
     char** flags = tune16_.flags;

     if ((strncmp(flags[0], "PGAUSS", 5) == 0) ||
         (strncmp(flags[0], "GAUSS",  5) == 0)) {
          sprintf(OUTSTRING,"%10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n",
          matrix[0], matrix[1], matrix[2],
          matrix[3], matrix[4], matrix[5], matrix[6]); 
     }
     if ((strncmp(flags[0], "PGALTMODEL", 5) == 0) ||
         (strncmp(flags[0], "GALTMODEL",  5) == 0)) {
          sprintf(OUTSTRING,"%10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n",
          matrix[0], matrix[1], matrix[2],
          matrix[3], matrix[4], matrix[5], matrix[6], 
          matrix[7], matrix[8], matrix[9], matrix[10]); 
     }
     if (strncmp(flags[0], "SERSIC", 5) == 0){ 
          sprintf(OUTSTRING,"%10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n",
          matrix[0], matrix[1], matrix[2],
          matrix[3], matrix[4], matrix[5], matrix[6], 
          matrix[7]); 
     }
     if (strncmp(flags[0], "EXTPGAUSS", 5) == 0){ 
          sprintf(OUTSTRING,"%10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n",
          matrix[0], matrix[1], matrix[2],
          matrix[3], matrix[4], matrix[5], matrix[6], 
          matrix[7], matrix[8]); 
     }

     return OUTSTRING;
}


