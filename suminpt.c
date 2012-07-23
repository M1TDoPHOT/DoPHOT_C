#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tuneable.h"
#include "suminpt.h"

/* dophot subroutine converted to c void function 01-25-2012 */

int suminpt_(int* I1, int* I2, float* STARPAR, int* NPAR, char* objectline)
{

     int nitems, noerr; 
     char** flags = tune16_.flags;

     if (*NPAR < 7){
          fprintf(stderr, "Less than 7 parameters allocated in SUMINPT routine... Dire warning !!!\n") ;
     }

     if ((strncmp(flags[0], "PGAUSS", 5) == 0) ||
         (strncmp(flags[0], "GAUSS",  5) == 0)) {
          nitems = sscanf(objectline, "%d %d %e %e %e %e %e %e %e",
                            I1, I2 , STARPAR+0, STARPAR+2,
                                     STARPAR+3, STARPAR+1,
                                     STARPAR+4, STARPAR+5,
                                     STARPAR+6);
          noerr = (nitems == 9);
     }
     if ((strncmp(flags[0], "PGALTMODEL", 5) == 0) ||
         (strncmp(flags[0], "GALTMODEL",  5) == 0)) {
          nitems = sscanf(objectline, "%d %d %e %e %e %e %e %e %e %e %e %e %e",
                            I1, I2 , STARPAR+0, STARPAR+2,
                                     STARPAR+3, STARPAR+1,
                                     STARPAR+4, STARPAR+5,
                                     STARPAR+6, STARPAR+7,
                                     STARPAR+8, STARPAR+9,
                                     STARPAR+10);
          noerr = (nitems == 13);
     }
     if (strncmp(flags[0], "SERSIC", 5) == 0){
          nitems = sscanf(objectline, "%d %d %e %e %e %e %e %e %e %f",
                            I1, I2 , STARPAR+0, STARPAR+2,
                                     STARPAR+3, STARPAR+1,
                                     STARPAR+4, STARPAR+5,
                                     STARPAR+6, STARPAR+7);
          noerr = (nitems == 10);
     }
     if (strncmp(flags[0], "EXTPGAUSS", 5) == 0){
          nitems = sscanf(objectline, "%d %d %e %e %e %e %e %e %e %f %f",
                            I1, I2 , STARPAR+0, STARPAR+2,
                                     STARPAR+3, STARPAR+1,
                                     STARPAR+4, STARPAR+5,
                                     STARPAR+6, STARPAR+7,
                                     STARPAR+8);
          noerr = (nitems == 11);
     }
 
     return noerr;
}
