#include <stdio.h>
#include <math.h>
#include "tuneable.h"
#include "estuff_struct.h"
#include "eoff_struct.h"
#include "string.h"
#include "cast_arr.h"
#include "ellipse.h"
#include "outputs.h"

/* dophot subroutine converted to c void function 02-06-2012 */

char* stdotpt_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPARR_ptr, float* APPARR, int* lNAPPLE_ptr, float* PROBG_ptr, int which_model)
{

     /* dereference the pointers */
     int   I      = *I_ptr;
     int   ITYPE  = *ITYPE_ptr;
     int   lNAPPLE = *lNAPPLE_ptr;
     float PROBG  = *PROBG_ptr; 
     
     /* renaming used common block variables */
     char** flags = tune16_.flags;
     float XFRAC = eoff_.xfrac;
     float YFRAC = eoff_.yfrac;
     short int* EMSUB = estuff_.emsub;
     float*     EMERR = estuff_.emerr;
     float**    EMPAR = estuff_.empar;

     /* substace of subroutine begins here */
     /* Get image shape and area and calc. Fit magnitude: */
     static char OUTSTRING[300];

     float AREA, AMAJOR, AMINOR, TILT;
     float FMAG;
     float APMAG, APERUNC;
     float XC, YC;
     float t1, t2, t3, t4;

     if (ITYPE != 8){
          ellipse_(STPARR+4, STPARR+5, STPARR+6, 
                    &AREA, &AMAJOR, &AMINOR, &TILT);
          FMAG = AREA*STPARR[1];
          if (FMAG <= 0.0f){
               FMAG =  99.999f;
          }
          else{
               FMAG = -2.5f*log10f(FMAG);
          }
          TILT = 57.29578f*TILT;
     }
     else{
          if (STPARR[5] == -1.0f){
               FMAG =  99.999f;
          }
          else{
               FMAG = -99.999f;
          }
          AMAJOR = STPARR[4];
          AMINOR = STPARR[6];
          TILT   = STPARR[5];
//          if (STPARR[4] >= STPARR[6]){
//               AMAJOR = STPARR[4];
//               AMINOR = STPARR[6];
//               TILT   = 0.0f;
//          }
//          else{
//               AMAJOR = STPARR[6];
//               AMINOR = STPARR[4];
//               TILT   = 90.0f;
//          }
     }

     /* Convert aperture flux to magnitudes: */
     if (APPARR[0] <= 0.0f){
          APMAG = 99.999f;
     }
     else{
          APMAG = -2.5f*log10f(APPARR[0]);
     } 

     /* If available, get uncertainty in aperture magnitudes: */
     APERUNC = 99.999f;
     if (lNAPPLE >= 5){
          APERUNC = APPARR[4];
     }

     /* Fix the co-ordinates so that the center of the first pixel is 
        0.5 and not 1.0 as in the internal representation */
     XC = STPARR[2] - 0.5f;
     YC = STPARR[3] - 0.5f;

     /* process the empirical fit data: */
     I = I-1;
     if ((int)(EMSUB[I]) != 0){
          t1 = EMPAR[I][0];
          t2 = EMPAR[I][1]*10000.0f*AREA;
          if (t2 <= 0.0f){
               t2 = 99.999f;
          }
          else{
               t2 = -2.5*log10f(t2);
          }
          t3 = EMPAR[I][2] + XFRAC - 0.5f;
          t4 = EMPAR[I][3] + YFRAC - 0.5f;
     }
     else{
          t1 = EMPAR[I][0];
          t2 = -99.999f;
          t3 = EMPAR[I][2];
          t4 = EMPAR[I][3];
     }

     if ((strncmp(flags[0], "PGAUSS", 5) == 0) ||
         (strncmp(flags[0], "GAUSS",  5) == 0)) {
          sprintf(OUTSTRING," %4d %2d %8.2f %8.2f %8.3f %6.3f %9.2f %9.3f %9.3f %7.2f %10.2E %8.3f %6.3f %9.2f %7.3f %9.2f %8.3f %8.2f %8.2f %6.3f   \n", 
          I+1, ITYPE, XC, YC, FMAG, 
          APPARR[3], STPARR[0], AMAJOR, AMINOR, TILT, PROBG,
          APMAG, APERUNC, APPARR[1], APPARR[2], t1, t2, t3, t4, EMERR[I]);
          /* Order is:  No., obtype, xpos, ypos, fitmag, 
          err_fitmag, fitsky, FWHM_major, FWHM_minor, Tilt, probgal,
          apmag, err_apmag, apsky, diff_fit_ap */
     }
     if ((strncmp(flags[0], "PGALTMODEL", 5) == 0) ||
         (strncmp(flags[0], "GALTMODEL",  5) == 0)) {
          sprintf(OUTSTRING," %4d %2d %8.2f %8.2f %8.3f %6.3f %9.2f %9.3f %9.3f %7.2f %10.3E %10.3E %10.3E %10.3E %2d %10.2E %8.3f %6.3f %9.2f %7.3f %9.2f %8.3f %8.2f %8.2f %6.3f   \n", 
          I+1, ITYPE, XC, YC, FMAG, 
          APPARR[3], STPARR[0], AMAJOR, AMINOR, TILT,
          STPARR[7], STPARR[8], STPARR[9], STPARR[10], which_model, PROBG,
          APMAG, APERUNC, APPARR[1], APPARR[2], t1, t2, t3, t4, EMERR[I]);
          /* Order is:  No., obtype, xpos, ypos, fitmag, 
          err_fitmag, fitsky, FWHM_major, FWHM_minor, Tilt, 
          NEWI1, NEWI2, NEWI3, NEWI4, which_model, probgal,
          apmag, err_apmag, apsky, diff_fit_ap */
     }
     if (strncmp(flags[0], "SERSIC", 5) == 0){
          sprintf(OUTSTRING," %4d %2d %8.2f %8.2f %8.3f %6.3f %9.2f %9.3f %9.3f %7.2f %6.2f %2d %10.2E %8.3f %6.3f %9.2f %7.3f %9.2f %8.3f %8.2f %8.2f %6.3f   \n", 
          I+1, ITYPE, XC, YC, FMAG, 
          APPARR[3], STPARR[0], AMAJOR, AMINOR, TILT,
          STPARR[7], which_model, PROBG,
          APMAG, APERUNC, APPARR[1], APPARR[2], t1, t2, t3, t4, EMERR[I]);
          /* Order is:  No., obtype, xpos, ypos, fitmag, 
          err_fitmag, fitsky, FWHM_major, FWHM_minor, Tilt, 
          sersic index, which_model, probgal,
          apmag, err_apmag, apsky, diff_fit_ap */
     }
     if (strncmp(flags[0], "EXTPGAUSS", 5) == 0){
          sprintf(OUTSTRING," %4d %2d %8.2f %8.2f %8.3f %6.3f %9.2f %9.3f %9.3f %7.2f %6.2f %6.2f %2d  %10.2E %8.3f %6.3f %9.2f %7.3f %9.2f %8.3f %8.2f %8.2f %6.3f   \n", 
          I+1, ITYPE, XC, YC, FMAG, 
          APPARR[3], STPARR[0], AMAJOR, AMINOR, TILT,
          STPARR[7], STPARR[8], which_model, PROBG,
          APMAG, APERUNC, APPARR[1], APPARR[2], t1, t2, t3, t4, EMERR[I]);
          /* Order is:  No., obtype, xpos, ypos, fitmag, 
          err_fitmag, fitsky, FWHM_major, FWHM_minor, Tilt, 
          b4, which_model, probgal,
          apmag, err_apmag, apsky, diff_fit_ap */
     }

     return OUTSTRING;
}


