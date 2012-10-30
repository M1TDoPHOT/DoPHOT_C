#include <stdio.h>
#include "free_parking_struct.h"
#include "cast_arr.h"
#include "errupd.h"

/* dophot subroutine converted to c void function on 01-17-2012 */

void errupd_(float *C_ptr, float *SHADERR, int *NFIT_ptr)
{
     /* dereferencing the pointers */
     int NFIT = *NFIT_ptr ;
     /* SHADDERR is a pointer to the first element of a 1D array and
        should not be de referenced       */
     /* recasting the 2 dimensional array */
     float** C = free_parking_.npmaxbynpmaxarray;
     recast_float_1dto2darr(NFIT, NFIT, C_ptr, C); 

     /* substance of the suboutine begins here */
     int I ;
     for (I = 0; I < NFIT; I++){
          SHADERR[I] = C[I][I]*C[I][I] ;
     }

     /* no recasting needed as the only changed variable is 
        cast in both fortran and c as a 1d arr */

} 
