#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tuneable.h"
#include "subraster_struct.h"
#include "crudestat_struct.h"
#include "fitarrays_struct.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "parinterp.h"
#include "chisq.h"
#include "twofit.h"

/* dophot function converted to c void functoin 02-22-2012 */
// fits two stars initially by fixing shape parameters and finding positions

//twostar is a pointer to the function name TWOSTAR which returns double
double twofit_( double (*TWOSTAR)(short int*, float*, float*, int*, int*), float* STARPAR)
{
     /* renaming used (and pointer passed) commonblock variables */
     short int* IRECT = tune2_.irect;
     float* AVA  = tune15_.ava; //average shape values
     float* ACC  = tune15_.acc; //accuracies for the fit
     float* ALIM = tune15_.alim; //limits for the fit
     float* B   = fitarrays_.b;
     float* A   = fitarrays_.a;
     //the following are passed to chisq
     int NFIT2 = tune4_.nfit2 ; //number of shape parameters
     int NIT2  = 2 * tune4_.nit; //number of iterations
     int NPT   = crudestat_.npt;
     short int** XX = subraster_.xx; //(x,y) positions
     float*  Z = subraster_.z ; //intensities
     float* YE = subraster_.ye; //errors on intensities
     float* FB = fitarrays_.fb;
     float* C_ptr = fitarrays_.c ; //not recast, only passed

     /* substance of function begins here */
     float* BACC = malloc_float_1darr(2*NPMAX);
     float* BLIM = malloc_float_1darr(2*NPMAX);
     int BADFIT = 0; //false
     int CONV   = 1; //true
     double parinterp_return;
     double chisq_return;
     //number of parameters to fit (sky + 3 for each star)
     int seven_dum = 7; 
     float A4, A6;
     float ANGLE, ROOT, ROOT1, ROOT2, PMET;
     float DX2, DY2, DX, DY;
     float GAREA, SAREA;
     float DX74, DY74, DXMAX, DYMAX;
     float DOT, FAC1, FAC2;
     int I;
     float temp;
     float TWOFIT;

     int model = 0; //default 0 for pseudogauss and gauss
     if (strncmp(tune16_.flags[0], "EXTPGAUSS", 5) == 0){
          model = 1;
     }

     /* SENSE IS THAT IF X IS LONG, ANGLE IS SMALL. */
     /* SENSE IS THAT IF X&Y ARE POSITIVELY CORRELATED, ANGLE IS POSITIVE. */

     /* Changed indices. */
     parinterp_return = parinterp_(STARPAR+2, STARPAR+3, AVA);
     if (B[0] == 0.0f){
          A4 = 1.0f/A[4];
          A6 = 1.0f/A[6];
          ANGLE = atan2f(-2.0f*A[5], A6 - A4)/2.0f;
          PMET  = (A4 - A6)*(A4 - A6) + 4.0f*A[5]*A[5];
          ROOT  = sqrtf(PMET);
          ROOT1 = (A4 + A6 + ROOT)/2.0f;
          ROOT2 = ROOT1 - ROOT;
          DX2   = A[4] - AVA[4];
          DY2   = A[6] - AVA[6];
          DX    = sqrtf(max1(DX2,0.0f))/2.0f;
          DY    = sqrtf(max1(DY2,0.0f))/2.0f;
          BADFIT = ( max1(DX, DY) == 0.0f );
          DX    = fabsf(DX) * cosf(ANGLE)/fabsf(cosf(ANGLE)) ;
          DY    = fabsf(DY) * sinf(ANGLE)/fabsf(sinf(ANGLE)) ;
          GAREA = 1.0f/sqrtf( fabsf(ROOT1*ROOT2)   );
          SAREA = sqrtf( fabsf(AVA[4]*AVA[6]) );

          /* Changed indices. */
          DX74 = STARPAR[2] - A[2];
          DY74 = STARPAR[3] - A[3];
         
          DOT = DX74*DX + DY74*DY; 
          if (DOT > 0.0f){
               FAC1 = 0.6666666f;
	       FAC2 = 1.3333333f;
          }
          else{
	       FAC1 = 1.3333333f;
	       FAC2 = 0.6666666f;
          }

          /* Changed indices! */
          /* ADDING 3 EXTRA INDICES FOR X,Y,INTENSITY OF SECOND OBJ */
          B[0] = A[0];
          B[1] = A[1]*GAREA/SAREA/2.0f*FAC2;
          B[2] = A[2] - DX*FAC1;
          B[3] = A[3] - DY*FAC1;
          B[4] = A[1]*GAREA/SAREA/2.0f*FAC1;
          B[5] = A[2] + DX*FAC2;
          B[6] = A[3] + DY*FAC2;

     }// end if B[0] == 0
     
     /* Changed indices. */
     B[1] = logf(B[1]); //image z's moved to log space for fitting
     B[4] = logf(B[4]); //image z's moved to log space for fitting
     B[7] = AVA[4];
     B[8] = AVA[5];
     B[9] = AVA[6];
     if (NFIT2 > 7){ //non pseudogaussian variables set to avg variables
          for(I = 7; I < NFIT2; I++){
               B[I+3] = AVA[I];
          }
     } 
     if (model == 1){
          B[10] = logf(B[10]);
          B[11] = logf(B[11]);
     }

     for(I = 0; I < 4; I++){
          BACC[I+3] = ACC[I];
          BACC[I]   = ACC[I];
          BLIM[I+3] = ALIM[I];
          BLIM[I]   = ALIM[I];
     }

     /* Changed indices. */
     BACC[1] = -0.01f; //logarathimic limits of intensities
     BACC[4] = -0.01f;
     BLIM[1] = -10.0f;
     BLIM[4] = -10.0f;
     DXMAX = max1(fabsf(B[2]), fabsf(B[5]));
     DYMAX = max1(fabsf(B[3]), fabsf(B[6]));

     BADFIT = (BADFIT || (DXMAX > (float)(IRECT[0])/2.0f));
     BADFIT = (BADFIT || (DYMAX > (float)(IRECT[1])/2.0f));
     if (!BADFIT){
          chisq_return = chisq_(TWOSTAR, XX, Z, YE, &NPT, B, FB, C_ptr, 
                          &seven_dum, BACC, BLIM, &NIT2);
          TWOFIT = (float)chisq_return;
          CONV   = (TWOFIT < 1.0e10f);
     }

     /* Changed indices. */
     DXMAX = max1(fabsf(B[2]), fabsf(B[5]));
     DYMAX = max1(fabsf(B[3]), fabsf(B[6]));
     BADFIT = (BADFIT || (DXMAX > (IRECT[0]/2.0f + 1.0f)) );
     BADFIT = (BADFIT || (DYMAX > (IRECT[1]/2.0f + 1.0f)) );
     if (CONV){
          /* Changed indices. */
          // making the brighter one first
          if (B[1] < B[4]){
               for (I = 1; I < 4; I++){
                    temp   =  B[I];
                    B[I]   =  B[I+3];
                    B[I+3] = temp;
               }
          }
          /* Changed indices. */
          B[1] = expf(B[1]);
          B[4] = expf(B[4]);
          if (model == 1){
               B[10] = expf(B[10]);
               B[11] = expf(B[11]);
          }


          //final output array must read Aa[0-NPAR] Ab[0-NPAR]
          //array currently reads A[0] Aa[123] Ab[123] A[4-NPAR]
         
          B[NFIT2  ] = B[0]; //sky
          B[NFIT2+1] = B[4]; //intensity
          B[NFIT2+2] = B[5]; //x
          B[NFIT2+3] = B[6]; //y
          for(I = 4; I < NFIT2; I++){ //again, sigma + non-pgauss variables set to ava
               B[I      ] = AVA[I];
               B[I+NFIT2] = AVA[I];
          }

     }//end if CONV

     if ( (!CONV) || BADFIT){
          TWOFIT = 1.0e20f;
     }

     /* repointing possibly change pointers and freeing allocated memory */
     crudestat_.npt = NPT;
     free(BACC);
     free(BLIM);

     return TWOFIT;
}
