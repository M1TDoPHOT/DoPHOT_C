#include <stdlib.h>
#include <stdio.h>
#include "string.h"
#include "math.h"
#include "ellipse.h"
#include "outputs.h"

/* dophot subroutine converted to c void function 01-26-2012 */

char* badotpt_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPAR_ptr, float* APPARR, int* NAPPLE_ptr, float* PROBG_ptr)
{
     /* dereference pointers */
     int I       = *I_ptr      ;
     int ITYPE   = *ITYPE_ptr  ;
//unused     int NPAR   = *NPAR_ptr  ;
     int NAPPLE  = *NAPPLE_ptr ;
     float PROBG = *PROBG_ptr  ;

     /* substance of subroutine begins here */
     static char OUTSTRING[200];

     float FMAG, APMAG, APERUNC, XC, YC ; 
     float AREA, AMAJOR, AMINOR, TILT   ;

     /* Get image shape and area and calc. Fit magnitude: */
     if (ITYPE != 8){
          /* fortran function passes by pointers so c does for
             compatibility */
          ellipse_((STPARR+4),(STPARR+5),(STPARR+6),
                   &AREA, &AMAJOR, &AMINOR, &TILT);

          FMAG = AREA*STPARR[1] ;
          if (FMAG <= 0.0f) {
               FMAG = 99.999f ;
          }
          else{
               FMAG = -2.5f*log10f(FMAG) ;
          }
          TILT = 57.29578f * TILT ;
     }
     else{
          if (STPARR[5] == -1.0f){
               FMAG =  99.999f ;
          }
          else{
               FMAG = -99.999f ;
          }
          if (STPARR[4] >= STPARR[6]){
               AMAJOR = STPARR[4] ;
               AMINOR = STPARR[6] ;
               TILT   = 0.0f      ;
          }
          else{
               AMAJOR = STPARR[6] ;
               AMINOR = STPARR[4] ;
               TILT   = 90.0f     ;
          }
     }
           
     /* Convert aperture flux to magnitudes: */
     if (APPARR[0] <= 0.0f){
          APMAG = 99.999f ;
     }
     else{
          APMAG = -2.5f * log10f( APPARR[0] ) ;
     }

     /* If available, get uncertainty in aperture magnitudes: */
     APERUNC = 99.999f ;
     if (NAPPLE >= 5){
          APERUNC = APPARR[4];
     }

     /* Convert co-ordinates so that internal representation where 
        center of 1st pixel is 1.0 is changed so that center of 1st
        pixel is 0.5 */
     XC = STPARR[2] - 0.5f ;
     YC = STPARR[3] - 0.5f ;

     sprintf(OUTSTRING, "%6d%9.2f%9.2f%9.3f%9.3f%9.3f%8d.%9.2f%9.3f  \n",
                         I,  XC,   YC,   
                         FMAG, APPARR[3], STPARR[0], 
                         ITYPE, PROBG, APPARR[2]) ;

     return OUTSTRING;
}
