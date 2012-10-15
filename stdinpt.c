#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tuneable.h"
#include "math.h"
#include "ellipint.h"
#include "stdinpt.h"

/* dophot subroutine converted to c void function 01-25-2012 */

int stdinpt_(int* I1, int* I2, float* STPARR, int* NPAR, char* INSTR)
{
 
     /* substance of subroutine begins here */
     char** flags = tune16_.flags;
     float XC, YC, FMAG, AREA, AMAJOR, AMINOR, TILT, dum;
     float SKY, NEWI1, NEWI2, NEWI3, NEWI4, SERSIC_INDX, B4, B6;
     int nread, noerr;
  
     if (*NPAR < 7){
          fprintf(stderr, "Less than 7 parameters allocated in STDINPT routine... Dire warning !!!\n") ;
     }

     if ((strncmp(flags[0], "PGAUSS", 5) == 0) ||
         (strncmp(flags[0], "GAUSS",  5) == 0)) {
          nread = sscanf(INSTR, "%d        %d        %f        %f   \
                                 %f        %f        %f             \
                                 %f        %f        %f            ",
                                 I1,       I2,       &XC,      &YC, 
                                 &FMAG,    &dum,     &SKY,
                                 &AMAJOR,  &AMINOR,  &TILT);     
     }
     if ((strncmp(flags[0], "PGALTMODEL", 5) == 0) ||
         (strncmp(flags[0], "GALTMODEL",  5) == 0)) {
          nread = sscanf(INSTR, "%d        %d        %f        %f   \
                                 %f        %f        %f             \
                                 %f        %f        %f             \
                                 %g        %g        %g        %g   ",
                                 I1,       I2,       &XC,      &YC, 
                                 &FMAG,    &dum,     &SKY,
                                 &AMAJOR,  &AMINOR,  &TILT,
                                 &NEWI1,   &NEWI2,   &NEWI3,   &NEWI4);    
          STPARR[7]  = NEWI1 ;
          STPARR[8]  = NEWI2 ;
          STPARR[9]  = NEWI3 ;
          STPARR[10] = NEWI4 ;
     }
     if (strncmp(flags[0], "SERSIC", 5) == 0){
          nread = sscanf(INSTR, "%d        %d        %f        %f   \
                                 %f        %f        %f             \
                                 %f        %f        %f             \
                                 %f                                 ",
                                 I1,      I2,        &XC,      &YC, 
                                 &FMAG,    &dum,     &SKY,
                                 &AMAJOR,  &AMINOR,  &TILT,
                                 &SERSIC_INDX);    
          STPARR[7] = SERSIC_INDX ;
     }
     if (strncmp(flags[0], "EXTPGAUSS", 5) == 0){
          nread = sscanf(INSTR, "%d        %d        %f        %f   \
                                 %f        %f        %f             \
                                 %f        %f        %f             \
                                 %f        %f                       ",
                                 I1,      I2,        &XC,      &YC, 
                                 &FMAG,    &dum,     &SKY,
                                 &AMAJOR,  &AMINOR,  &TILT,
                                 &B4,      &B6);    
          STPARR[7] = B4 ;
          STPARR[8] = B6 ;
     }
     STPARR[0] = SKY ;
     noerr = (nread > 0);

     /* Convert position centers so that the center of the 1st pixel is 1.0 
        ( as required in the internal representation) rather than
        0.5 as used in the COMPLETE style output */

     STPARR[2] = XC + 0.5f ;
     STPARR[3] = YC + 0.5f ;

     if (*I2 == 8){
          STPARR[1] = 0.0f ;   
          STPARR[4] = AMAJOR ;
          STPARR[5] = TILT ;   
          STPARR[6] = AMINOR ;
//          if ((TILT > 89.5f) && (TILT < 90.5f)){
//               STPARR[1] = 0.0f ;  
//               STPARR[4] = AMINOR ;   
//               STPARR[5] = 0.0f ;   
//               STPARR[6] = AMAJOR ;   
//          }
//          else{
//               STPARR[1] = 0.0f ;   
//               STPARR[4] = AMAJOR ;
//               STPARR[5] = TILT ;   
//               STPARR[6] = AMINOR ;
//          }
     }
     else{
          /* fortran passes all values by pointers, so c versions do the same
             for compatibility */
          ellipint_(&AMAJOR,&AMINOR,&TILT,&AREA,STPARR+4,STPARR+5,STPARR+6) ;
          if (FMAG > 0.0f){
               STPARR[1] = 0.0f ;
          }
          else{
               STPARR[1] = powf(10.0f,-0.4f*FMAG) / AREA ;
          }
     }
     
     return noerr;
}
