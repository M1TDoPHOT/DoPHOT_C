#include <stdlib.h>
#include <math.h>
#include "tuneable.h"
#include "drfake_struct.h"
#include "mini_mathlib.h"
#include "cast_arr.h"
#include "gauss.h"

/* dophot function with 'entry' converted to 2 c double functions
   02-17-2012 */

double gauss2d_(short int* IX, float* A, float* FA, int* M_ptr, int* fitcall_ptr)
{

     // A[0] = SKY,
     // A[1] = MAX intensity
     // A[2], [3], x, y, positions in field
     // A[4], [5], [6], sigmax^2, sigmaxy, sigmay^2

     // if call from chisq, the parameters may be in a different
     // space, i.e. log space and will need to be converted.
     // else, the parameters will be in linear space.
     // fitcall == 1 means gaussian2d was called from chisq
     int fitcall = *fitcall_ptr;

     /* rename used common block variables */
     int   NEEDIT = drfake_.needit;

     /* substance of function begins here */
     double X,Y;
     double sigx, sigy;
     double Z, A1;
     double DIDZ, PEXP;
     double GAUSSIAN2D;

     double HALF   = 0.5;
     double EXPMIN = -23.0;

     /* Changed indices */
     X = (double)(IX[0]) - A[2];
     Y = (double)(IX[1]) - A[3];

     sigx = 1.0/A[4];
     sigy = 1.0/A[6];
     Z = HALF*(sigx*X*X + 2.0*A[5]*Y*X + sigy*Y*Y); 

     if (Z < 0.0){
          Z    = max1(Z, EXPMIN);
     }
     PEXP = exp(-Z);

     /* Changed indices */
     if (fitcall == 1){
          A1    = exp(A[1]); //log case
          FA[1] = A1*PEXP; //log case
     }
     else{
          A1    = A[1]; //linear case
          FA[1] = PEXP; //linear case
     }
     DIDZ = A1*PEXP;
     GAUSSIAN2D = A1*PEXP + A[0];

     if (NEEDIT){
          /* Changed indices */
          FA[0] = 1.0                       ;
          FA[2] = DIDZ*(sigx*X + A[5]*Y)    ;
          FA[3] = DIDZ*(A[5]*X + sigy*Y)    ;
          FA[4] = DIDZ*(HALF*sigx*sigx*X*X);
          FA[5] = DIDZ*(-X*Y)               ;
          FA[6] = DIDZ*(HALF*sigy*sigy*Y*Y);
     }

     return GAUSSIAN2D;
}

double gauss4d_(short int* IX, float* A, float* FA, int* M_ptr, int* fitcall_ptr)
{

     // A[0] = SKY,
     // A[1] = MAX intensity for obj 1
     // A[2], [3], x, y, positions in field for obj 1
     // A[4] = MAX intensity for obj 1
     // A[5], [6], x, y, positions in field for obj 1
     // A[7], [8], [9], sigmax^2, sigmaxy, sigmay^2 for both gaussians
     
     /* substance of function begins here */
     double* PP = malloc(2*sizeof(double));
     PP[0] = 0.0;
     PP[1] = 0.0;
     int   I, IOFF;
     double X,Y;
     double sigx, sigy;
     double Z, A1;
     double DIDZ;
     double GAUSSIAN4D;

     double HALF   = 0.5;

     sigx = 1.0/A[7];
     sigy = 1.0/A[9];    
     for (I = 0; I < 2; I++){
          IOFF = 3*I;
          /* changed indices */
          X = (double)(IX[0]) - A[2 + IOFF];
          Y = (double)(IX[1]) - A[3 + IOFF];
          Z = HALF*(sigx*X*X + 2.0*A[8]*Y*X + sigy*Y*Y); 
          PP[I] = exp(-Z);

          // default caled by chisq
          A1           = exp(A[1 + IOFF]);
          FA[1 + IOFF] = A1*PP[I];
          
          PP[I] = A1*PP[I];
          DIDZ  = A1*PP[I];

          /* Changed indices and A3 to A1. */
          FA[2 + IOFF] = DIDZ*(sigx*X + A[8]*Y);
          FA[3 + IOFF] = DIDZ*(A[8]*X + sigy*Y);
     }

     FA[0]   = 1.0;
     GAUSSIAN4D = PP[0] + PP[1] + A[0]; 

     free(PP);
     return GAUSSIAN4D;
}


