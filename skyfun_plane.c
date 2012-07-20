#include "undergnd_struct.h"
#include "skyfun_plane.h"
/* dophot function converted to C function 01-20-2012 */

double skyfun_(short int* IX, float* A, float* FA, int* M_ptr, int* MMAX_ptr)
{
 
     /* M and MMAX are never used and so never need be dereferecned */
     /* the rest are 1d arrays and are dereferenced by matrix notation */
     float X = (float)IX[0];
     float Y = (float)IX[1];

     /* renaming the global variables */
     int NFAST  = undergnd_.NFAST ;
     int NSLOW  = undergnd_.NSLOW ;

     /* substance of subroutine starts here */
     FA[0] = 1.0;
     FA[1] = 0.5*(X - 0.5*NFAST)/(0.5*NFAST) ;
     FA[2] = 0.5*(Y - 0.5*NSLOW)/(0.5*NSLOW) ;

     float SKYF = 0.0 ;
     int I            ;
     for (I = 0; I < 3; I++){
          SKYF = SKYF + A[I]*FA[I] ;
     }
     return SKYF ;

}
