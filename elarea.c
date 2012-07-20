#include <math.h>
#include "logh.h"
#include "tuneable.h" 
#include "elarea.h"

/* dophot function converted to c function 01-29-2012 */

double elarea_(float* B5_ptr, float* A6_ptr, float* B7_ptr)
{

     /* dereference the pointers */
     float B5 = *B5_ptr;
     float A6 = *A6_ptr;
     float B7 = *B7_ptr;
     /* renaming the used tuneable parameter */
     int lverb = tune14_.lverb;

     /* substance of subroutine begins here */
     float A5, A7, ROOT, ROOT1, ROOT2;
     double ret_val;

     A5    = 1.0f/B5;
     A7    = 1.0f/B7;
     ROOT  = sqrtf( (A5 - A7)*(A5 - A7) + (4.0f*A6*A6) );
     ROOT1 = (A5 + A7 + ROOT)/2.0f;
     ROOT2 = ROOT1 - ROOT;

     if (ROOT2 > 0.0f){
          ret_val = 1.0/sqrt(ROOT1*ROOT2);
     }
     else{
          if (lverb > 20){
               fprintf(logfile,"ROOTS 1 & 2 = %f %f (....ELAREA)\n", ROOT1, ROOT2);
          }
          ret_val = 0.0;
     }
     
     return ret_val;
}

