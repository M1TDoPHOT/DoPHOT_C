#include "drfake_struct.h"
#include "math.h"
#include "mini_mathlib.h"
#include "skyfun_hub.h"

/* dophot function converted to c double function 02-03-2012 */
/*
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:  NOTE THAT NONE OF THE QUANTITIES ARE SCALED.  IN PRINCIPAL WE COULD
C:  SCALE BOTH THE INTENSITIES AND THE POSITIONS, THOUGH THE LATTER LOOKS
C:  TRICKY TO ME.  NOTE THAT DERIVATIVES OF 5, 6 AND 7 ARE COMPUTED BUT
C:  NOT FIT FOR AT THE MOMENT.  WE ADOPT THE SHAPE OF THE SEEING PROFILE
C:  BUT COULD EASILY SOLVE FOR ALL 7 PARAMETERS IN VARIPAR.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/

double hubfun_(short int* IX, float* A, float* FA, int* M_ptr, int* MMAX_ptr)
{
    
     /* never used M_ptr or MMAX_ptr, so no need to dereference */
     /* rename used common block variables */
     int NEEDIT = drfake_.needit;
     
     /* substance of function begins here */
     float HALF   = 0.5f;
     float EXPMIN = -23.0f;

     float X,Y;
     float A5, A7;
     float T5, T6, T7, T1;
     float DENOM, DDDT, PEXP;
     float HUBFUN;

     X  = (float)IX[0] - A[1];
     Y  = (float)IX[1] - A[2];
     A5 = 1.0f/A[4];
     A7 = 1.0f/A[6];
     T5 = A5*X     ;
     T6 = A[5]*Y   ;
     T7 = A7*Y     ;
     T1 = HALF*((T5+2.0f*T6)*X + T7*Y);
     if (T1 > 0.0f){
          DENOM = 1.0f + T1  ;
          DDDT  = 1.0f       ;
          PEXP  = 1.0f/DENOM ;
     }
     else{
          T1    = max1(T1, EXPMIN);
          PEXP  = expf(-T1)        ;
          DENOM = 1.0f            ;
          DDDT  = 1.0f            ;
     }
     FA[3] = PEXP        ;
     PEXP = A[3]*PEXP    ;
     HUBFUN = PEXP + A[0];
     if (NEEDIT){
          FA[5] = PEXP*DDDT/DENOM    ;
          FA[1] = (T5 + T6)*FA[5]    ;
          FA[2] = (A[5]*X + T7)*FA[5];
          FA[4] = HALF*T5*T5*FA[5]   ;
          FA[6] = HALF*T7*T7*FA[5]   ;
          FA[5] = -X*Y*FA[5]         ;
          FA[0] = 1.0f               ;
     }

     return HUBFUN;
}
         
/* hub2fun, hubfun's second entry, was never called by dophot and thus not translated to c */ 

