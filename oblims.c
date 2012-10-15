#include <math.h>
#include <string.h>
#include "logh.h"
#include "tuneable.h"
#include "mini_mathlib.h"
#include "oblims.h"

/* dophot subroutine converted to c void function 02-02-2012 */
/* set JRECT size to reflect size of star ellipse (sigma[xy]) */
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   PSEUDOGAUSSIAN EXP(-T**2) = 1/(1 + T**2 + T**4/2 + T**6/6)
  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

void oblims_(float* STARPAR, short int* JRECT)
{

     /* both passed vars are 1d arrays and do not need to be dereferenced 
        as a seperate step */
     /* renaming tuneable variables */
     float NPHOB = (float)tune6_.nphob; 
     
     /* substance of subroutine begins here */
     float TEMP, FUDGEX, FUDGEY;
     int I;
     float BETA6;

     if (strncmp(tune16_.flags[0], "EXTPGAUSS", 5) == 0){
          BETA6 = STARPAR[8];
     }
     else{
          BETA6 = tune17_.beta6;
     }

     if (STARPAR[1] > 0.0f){
          fprintf(logfile, "starpar[1] = %9.6f\n", STARPAR[1]);
          TEMP = powf( ((STARPAR[1]/NPHOB)*(6.0f/BETA6)), 0.33333333f );

          if (STARPAR[4] > 0){
               FUDGEX = sqrtf(TEMP*STARPAR[4]*2.0f);
               fprintf(logfile, "fudgex, temp & starpar[4] = \n");
               fprintf(logfile, "%9.6f  %9.6f  %9.6f  \n", FUDGEX, TEMP, STARPAR[4]);
          }
          else{
               FUDGEX = 10.0f;
          }

          if (STARPAR[6] > 0.0f){
               FUDGEY = sqrtf(TEMP*STARPAR[6]*2.0f);
          }
          else{
               FUDGEY = 10.0f;
          }
     }
     else{
          FUDGEX = 10.0f;          
          FUDGEY = 10.0f;
     }
     JRECT[0] = STARPAR[2] - FUDGEX;
     JRECT[1] = STARPAR[2] + FUDGEX;
     JRECT[2] = STARPAR[3] - FUDGEY;
     JRECT[3] = STARPAR[3] + FUDGEY;

     fprintf(logfile,"JRECT[I], I = ");
     for (I = 0; I < 4; I++){
          fprintf(logfile,"%9d %9d  ", JRECT[I], I);
     }
     fprintf(logfile,"\n");

}
