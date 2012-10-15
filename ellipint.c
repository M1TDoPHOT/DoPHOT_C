#include <math.h>
#include "logh.h"
/* dophot subroutine converted to c void function 01-13-2012 */
/* converts from amajor, minor, tilt to sigmax^2, sigmay^2, sigmaxy */

void ellipint_(float *AMAJOR_ptr, float *AMINOR_ptr, float *TILT_ptr, float *AREA_ptr, float *B5_ptr, float *A6_ptr, float *B7_ptr)
{
     /* dereferencing the pointers */
     float AMAJOR  = *AMAJOR_ptr ;
     float AMINOR  = *AMINOR_ptr ;
     float TILT    = *TILT_ptr   ;
     float AREA    = *AREA_ptr   ;
     float B5      = *B5_ptr     ;
     float A6      = *A6_ptr     ;
     float B7      = *B7_ptr     ;

     /* actual meat of the subroutine begins here */
     if (AMAJOR < AMINOR) {
          fprintf(logfile, " Maj-axis .LT. Min-axis in WARMSTART input data \n");
          fprintf(logfile, " Proceed at grave risk !!! \n");
     }
     if (AMAJOR <= 0.01f){
          AMAJOR = 0.01f ;
     }
     if (AMINOR <= 0.01f){
          AMINOR = 0.01f ;
     }

     float ROOT1, ROOT2;
     ROOT1 = (2.35482f/AMINOR)*(2.35482f/AMINOR) ;
     ROOT2 = (2.35482f/AMAJOR)*(2.35482f/AMAJOR) ;
     AREA  =  6.2832f / sqrtf(ROOT1*ROOT2)       ;
          
     float A7MNA5, A7PLA5, A5, A7  ;
     A6 = 0.5f*(ROOT2 - ROOT1) * sinf( 2.0f*(TILT/57.29578f) );
     A7MNA5 = (ROOT1 - ROOT2) * cosf( 2.0f*(TILT/57.29578f) );
     A7PLA5 = (ROOT1 + ROOT2)   ;
     A7 = 0.5f*(A7PLA5 + A7MNA5) ;
     A5 = 0.5f*(A7PLA5 - A7MNA5) ;

     B5 = 1.0f/A5 ;
     B7 = 1.0f/A7 ;

     /* rereferencing the pointers */
     *AMAJOR_ptr = AMAJOR ;
     *AMINOR_ptr = AMINOR ;
     *TILT_ptr   = TILT   ;
     *AREA_ptr   = AREA   ;
     *B5_ptr     = B5     ;
     *A6_ptr     = A6     ;
     *B7_ptr     = B7     ;
}
