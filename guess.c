#include "tuneable.h"
#include "crudestat_struct.h"
#include "unitize_struct.h"
#include "parinterp.h"

/* dophot function converted to c double function 02-05-2012 */
/*
C:  SERIOUS PROBLEM HERE.  I'D LIKE TO STOP MAKING GUESS1 A FUNCTION.
C:  BUT IT'S CALLED BY ISEARCH, FILLERUP AND DOPHOT.  SO I NEED TO MAKE
C:  CHANGES CONSISTENTLY
*/
/* guess 1 and guess 2 can also be combined, but with a toggle
   they only differ by one line */

double guess1_(float* A, float* STARPAR, int* IX_ptr, int* IY_ptr)
{

     /* dereference pointers */
     int IX = *IX_ptr; 
     int IY = *IY_ptr; 

     /* rename used common block vars */
     float XM      = crudestat_.xm;
     float YM      = crudestat_.ym;
     float SUM2    = crudestat_.sum2;
     float MAXVAL  = (float)crudestat_.maxval;
     float UFACTOR = unitize_.ufactor;

     /* substance of function begins here */
     float XT = (float)IX;
     float YT = (float)IY;
     float GUESS1;
     double dummy;
   
     dummy = parinterp_(&XT, &YT, A);
     GUESS1 = A[0];
     A[0] = SUM2/UFACTOR;
     A[1] = MAXVAL/UFACTOR - A[0];
     A[2] = XM;
     A[3] = YM;
  
     return GUESS1;
}

double guess2_(float* A, float* STARPAR, int* IX_ptr, int* IY_ptr)
{

     /* dereference pointers */
     int IX = *IX_ptr; 
     int IY = *IY_ptr; 

     /* rename used common block vars */
     int NPAR      = tune4_.npar;
     float UFACTOR = unitize_.ufactor;

     /* substance of function begins here */
     float GUESS2;
     int K;

     IX = (int)(STARPAR[2] + 0.5f);
     IY = (int)(STARPAR[3] + 0.5f);
     for (K = 4; K < NPAR; K++){
          A[K] = STARPAR[K];
     }
     A[0] = STARPAR[0]/UFACTOR;
     A[1] = STARPAR[1]/UFACTOR;
     A[2] = STARPAR[2] - IX;
     A[3] = STARPAR[3] - IY;
     GUESS2 = STARPAR[0];

     /* reference changed pointers */
     *IX_ptr = IX;
     *IY_ptr = IY;

     return GUESS2;
}

double guess3_(float* A, float* STARPAR, int* IX_ptr, int* IY_ptr)
{

     /* dereference pointers */
     int IX = *IX_ptr; 
     int IY = *IY_ptr; 

     /* rename used common block vars */
     float UFACTOR = unitize_.ufactor;

     /* substance of function begins here */
     float GUESS3;
     double dummy;

     IX = (int)(STARPAR[2] + 0.5f);
     IY = (int)(STARPAR[3] + 0.5f);
     dummy = parinterp_(STARPAR+2, STARPAR+3, A);
     A[0] = STARPAR[0]/UFACTOR;
     A[1] = STARPAR[1]/UFACTOR;
     A[2] = STARPAR[2] - IX;
     A[3] = STARPAR[3] - IY;
     GUESS3 = STARPAR[0];

     /* reference changed pointers */
     *IX_ptr = IX;
     *IY_ptr = IY;

     return GUESS3;
}
     
