#include "tuneable.h"
#include "unitize_struct.h"
#include "parupd.h"

/* dophot subroutine converted to c void function 02-05-2012 */
/*
C:  THIS SUBROUTINE IS PROFILE SPECIFIC AND NEEDS REWRITING FOR EACH
C:  FITTING FUNCTION IF A[0], A[1], A[2] or A[3] are changed.
C:  ELSE IT's FINE.
*/

void parupd_(float* A, float* STARPAR, int IX, int IY, int NPAR)
{

     /* rename used common block vars */
     float UFACTOR = unitize_.ufactor;

     /* substance of subroutine begins here */
     int K;
    
     for (K = 0; K < NPAR; K++){
          STARPAR[K] = A[K];
     }
     STARPAR[0] = A[0]*UFACTOR;
     STARPAR[1] = A[1]*UFACTOR;
     STARPAR[2] += (float)IX;
     STARPAR[3] += (float)IY;
    
} 

//twoupd currently unused
void twoupd_(float* A, float* STARPAR, int* IX_ptr, int* IY_ptr)
{

     /* dereference pointers */
     int IX = *IX_ptr;
     int IY = *IY_ptr;

     /* substance of subroutine begins here */
     int K;

     STARPAR[0] = A[0];
     for (K = 0; K < 4; K += 3){
          STARPAR[K + 1] = A[K + 1];
          STARPAR[K + 2] = A[K + 2] + (float)IX;
          STARPAR[K + 3] = A[K + 3] + (float)IY;
     }
    
} 
