#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tuneable.h"
#include "drfake_struct.h"
#include "mini_mathlib.h"
#include "cast_arr.h"
#include "pgauss.h"

/* dophot function with 'entry' converted to 2 c double functions
   02-17-2012 */

double pgauss2d_(short int* ix, float* a, float* fa, int* m_ptr, int* fitcall_ptr)
{
     // a[0] = sky,
     // a[1] = max intensity
     // a[2], [3], x, y, positions in field
     // a[4], [5], [6], sigmax^2, sigmaxy, sigmay^2

     // if call from chisq, the parameters may be in a different
     // space, i.e. log space and will need to be converted.
     // else, the parameters will be in linear space.
     // fitcall == 1 means pseud2d was called from chisq
     int fitcall = *fitcall_ptr;

     /* rename used common block variables */
     double beta4  = (double)tune17_.beta4;
     double beta6  = (double)tune17_.beta6;
     int   needit = drfake_.needit;

     /* substance of function begins here */
     double x,y;
     double sigx, sigy;
     double tx, ty, z, z2, a1;
     double denom, dddt, pexp;
     double pseud2d;

     double half   = 0.5;
     double third  = 0.3333333;
     double expmin = -23.0;

     /* changed indices */
     x = (double)(ix[0]) - a[2];
     y = (double)(ix[1]) - a[3];

     sigx = 1.0/a[4];
     sigy = 1.0/a[6];
     tx = sigx*x     ;
     ty = sigy*y     ;
     z = half*(tx*x + 2.0*a[5]*y*x + ty*y);

     if (z >= 0.0){
          z2    = half*z*z;
          denom = 1.0 + z + beta4*z2 + (beta6*third)*z*z2;
          dddt  = 1.0 + beta4*z + beta6*z2;
          pexp  = 1.0/denom;
     }
     else{
          z    = max1(z, expmin);
          pexp  = exp(-z);
          denom = 1.0;
          dddt  = 1.0;
     }
     
     /* changed indices */
     if (fitcall == 1){
          a1    = exp(a[1]); //log case
          fa[1] = a1*pexp; //log case
     }
     else{
          a1    = a[1]; //linear case
          fa[1] = pexp; //linear case
     }
     pexp  = a1*pexp;
     pseud2d = pexp + a[0];

     if (needit){
          /* changed indices */
          fa[5] = pexp*dddt/denom    ;
          fa[2] = (tx + a[5]*y)*fa[5]    ;
          fa[3] = (a[5]*x + ty)*fa[5];
          fa[4] = half*tx*tx*fa[5]   ;
          fa[6] = half*ty*ty*fa[5]   ;
          fa[5] = -x*y*fa[5]         ;
          fa[0] = 1.0                ;
     }
     return pseud2d;
}


double pgauss4d_(short int* ix, float* a, float* fa, int* m_ptr, int* lmmax_ptr)
{

     // a[0] = sky,
     // a[1] = max intensity for obj 1
     // a[2], [3], x, y, positions in field for obj 1
     // a[4] = max intensity for obj 1
     // a[5], [6], x, y, positions in field for obj 1
     // a[7], [8], [9], sigmax^2, sigmaxy, sigmay^2 for both gaussians
     

     /* rename used common block variables */
     double beta4  = (double)tune17_.beta4;
     double beta6  = (double)tune17_.beta6;

     /* substance of function begins here */
     double pp[2];
     pp[0] = pp[1] = 0;
     int   i, ioff;
     double x,y;
     double sigx, sigy, a1;
     double z, z2;
     double denom, dddt, fac;
     double pseud4d;

     double half   = 0.5;
     double third  = 0.3333333;

     sigx = 1.0/a[7];
     sigy = 1.0/a[9];    
     for (i = 0; i < 2; i++){
          ioff = 3*i;
          /* changed indices */
          x = (double)(ix[0]) - a[2 + ioff];
          y = (double)(ix[1]) - a[3 + ioff];
          z  = half*(sigx*x*x + 2.0*a[8]*y*x + sigy*y*y);
          if (z >= 0){
               z2    = half*z*z;
               denom = 1.0 + z + beta4*z2 + (third*beta6)*z*z2;
               dddt  = 1.0 + beta4*z + beta6*z2;
               pp[i] = 1.0/denom;
          }
          else{
               pp[i] = exp(-z);
               denom = 1.0f;
               dddt  = 1.0f;
          }

          /* changed indices and a3 to a1. */
          a1         = exp(a[1 + ioff]);
          fa[1 + ioff] = a1*pp[i];
          pp[i]        = a1*pp[i];

          fac = pp[i]*dddt/denom;

          /* changed indices */
          fa[2 + ioff] = (sigx*x + a[8]*y)*fac;
          fa[3 + ioff] = (a[8]*x + sigy*y)*fac;
     }
     fa[0]   = 1.0;
     pseud4d = pp[0] + pp[1] + a[0]; 

     return pseud4d;
}


