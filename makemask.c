#include <stdlib.h>
#include "tuneable.h"
#include "free_parking_struct.h"
#include "passmask_struct.h"
#include "passweird_struct.h"
#include "parinterp.h"
#include "cast_arr.h"
#include "makemask.h"

/* dophot subroutine converted to c void functoin 02-17-2012 */

//onestar is a pointer to the function name ONESTAR which returns double
//onestar is usually pseud2d
//void makemask_( double (*ONESTAR)(short int* IX, float* A, float* FA, int* M_ptr, int* fitcall_ptr) )
void makemask_( double (*ONESTAR)(short int*, float*, float*, int*, int*) )
{

     /* renaming used common block variables */
     int IXBY2 = tune12_.ixby2;
     int IYBY2 = tune12_.iyby2;
     int WEIRD = passweird_.weird;
     float** STARMASK = passmask_.starmask;

     /* substance of subroutine begins here */
     float* A  = free_parking_.npmaxarray_1;
     float* FA = free_parking_.npmaxarray_2;
     short int* IX = free_parking_.sifourarray_1;
     static float DUMX = 0.0f;
     static float DUMY = 0.0f;
     double parinterp_return;
     int I, J;
     int M;
     int ZERO_dum = 0;
     
     parinterp_return = parinterp_(&DUMX, &DUMY, A);

     /* Rearranged parameter values. */
     A[0] = 0.0f;
     A[1] = 1.0f;
     A[2] = 0.0f;
     A[3] = 0.0f;
     for (J = -IYBY2; J <= IYBY2; J++){
          IX[1] = J;
          for (I = -IXBY2; I <= IXBY2; I++){
               IX[0] = I;
               WEIRD = ( (I == 0) && (J == 0) );
               STARMASK[J + IYBY2][I + IXBY2] = (*ONESTAR)(IX, A, FA, &M, &ZERO_dum);
         }
     }

     /* reassigning renamed common block variables and freeing allocated mem */
     passweird_.weird =  WEIRD;
}
