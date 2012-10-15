#include <stdio.h>
#include <math.h>
#include "ellipse.h" 

/* dophot subroutine converted to c void function 01-13-2012 */
/* convert from sigmax^2, sigmay^2, sigmaxy to semimajor, semiminor, tilt, and area */

void ellipse_(float *B5_ptr, float *A6_ptr, float *B7_ptr, float *AREA_ptr, float *AMAJ_ptr, float *AMINO_ptr, float *ANGLE_ptr)
{
     /* dereferencing the pointers */
     float B5    = *B5_ptr    ;
     float A6    = *A6_ptr    ;
     float B7    = *B7_ptr    ;
     float AREA  = *AREA_ptr  ;
     float AMAJ  = *AMAJ_ptr  ;
     float AMINO = *AMINO_ptr ;
     float ANGLE = *ANGLE_ptr ;

     /* actual meat of the subroutine begins here */
     float A5 = 1.0f/B5 ;
     float A7 = 1.0f/B7 ;
     ANGLE = atan2f(-2.0f*A6, A7 - A5)/2.0f ;

     float ROOT, ROOT1, ROOT2 ;
     ROOT  = sqrtf( (A5 - A7)*(A5 - A7) + 4.0f*(A6*A6) ) ;
     ROOT1 = (A5 + A7 + ROOT) / 2.0f                     ;
     ROOT2 = ROOT1 - ROOT                                ;
     if (ROOT2 > 0.0f){
          AREA = 6.2832f/sqrt(ROOT1*ROOT2);
     }
     else{
	  AREA = 0                       ;
     }
     AMAJ  = 2.35482f*sqrtf(1.0f/ROOT2) ;
     AMINO = 2.35482f*sqrtf(1.0f/ROOT1) ;

     /* rereferencing the pointers */
     *B5_ptr = B5   ;
     *A6_ptr = A6   ;
     *B7_ptr = B7   ;
     *AREA_ptr = AREA   ;
     *AMAJ_ptr  = AMAJ  ;
     *AMINO_ptr = AMINO ;
     *ANGLE_ptr = ANGLE ;

}
