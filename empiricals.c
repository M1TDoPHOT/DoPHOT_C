#include <stdlib.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "drfake_struct.h"
#include "oimdev_struct.h"
#include "opqfake_struct.h"
#include "oefake_struct.h"
#include "orsfake_struct.h"
#include "odfake_struct.h"
#include "imdev_struct.h"
#include "pqfake_struct.h"
#include "efake_struct.h"
#include "rsfake_struct.h"
#include "dfake_struct.h"
#include "funnypass_struct.h"
#include "cast_arr.h"

/* dophot function converted to c function 02-07-2012 */

/* twelve point interpolation on both function and x and y derivatives */
/* a(3) and a(4) are the offsets of template (i.e. model) from center
   of current subraster. */
/* xy contains offset of current pixel from center of current raster */
/* note that since a(3) and a(4) change only rarely we can skip a
   lot of stuff. */
/* need to check boundaries of ix0 iy1 et cetera */
/* note that derivatives are now offset by 1/2 */

double empiricalinterp_(short int* xy, float* a, float* fa, int* m_ptr, int* fitcall_ptr, int whichemp)
{

     /* renaming used common block variables */
     /* which common block variables are used depends on
        whether oldemp or one emp was called */
     int NEEDIT = drfake_.needit;
     int ip, iq;
     float p, q;
     float e10, e01, e23, e32, e31, e13, e02, e20, e11, e22, e12, e21;
     int ir, is;
     float r, s;
     float d00, d01, d10, d11;
     float** dxmp;
     float** dymp;
     int** mp;
     int funny;


     /* for oldemp */
     if (whichemp == 0){
          ip = opqfake_.ip;
          iq = opqfake_.iq;
          p  = opqfake_.p;
          q  = opqfake_.q;

          e10 = oefake_.e10;
          e01 = oefake_.e01;
          e23 = oefake_.e23;
          e32 = oefake_.e32;
          e31 = oefake_.e31;
          e13 = oefake_.e13;
          e02 = oefake_.e02;
          e20 = oefake_.e20;
          e11 = oefake_.e11;
          e22 = oefake_.e22;
          e12 = oefake_.e12;
          e21 = oefake_.e21;

          ir = orsfake_.ir;
          is = orsfake_.is;
          r  = orsfake_.r;
          s  = orsfake_.s;

          d00 = odfake_.d00;
          d01 = odfake_.d01;
          d10 = odfake_.d10;
          d11 = odfake_.d11;
      
          dxmp = oimdev_.dxomp;
          dymp = oimdev_.dyomp;
          mp   = oimdev_.omp;

          funny = 0; // default false for oldemp;
     }
 
     /* for oneemp */
     if (whichemp == 1){
          ip = pqfake_.ip;
          iq = pqfake_.iq;
          p  = pqfake_.p;
          q  = pqfake_.q;

          e10 = efake_.e10;
          e01 = efake_.e01;
          e23 = efake_.e23;
          e32 = efake_.e32;
          e31 = efake_.e31;
          e13 = efake_.e13;
          e02 = efake_.e02;
          e20 = efake_.e20;
          e11 = efake_.e11;
          e22 = efake_.e22;
          e12 = efake_.e12;
          e21 = efake_.e21;

          ir = rsfake_.ir;
          is = rsfake_.is;
          r  = rsfake_.r;
          s  = rsfake_.s;

          d00 = dfake_.d00;
          d01 = dfake_.d01;
          d10 = dfake_.d10;
          d11 = dfake_.d11;
      
          dxmp = imdev_.dxemp;
          dymp = imdev_.dyemp;
          mp   = imdev_.emp;

          funny = funnypass_.funny;
     } 

     /* substance of function begins here */
     float flt_ip, flt_iq;
     float c00, c10, c01, c11; //four point coefficients
     float q0, p0, q1, p1; // handy definitions for symmetry
     int ix0, ix1, ix2, ix3;
     int iy0, iy1, iy2, iy3;
     int iu0, iv0, iu1, iv1;
     static int first = 1; //true
     static int second = 0; //false
     static float a2_oldemp = 0.0f;
     static float a3_oldemp = 0.0f;
     static float a2_oneemp = 1e10f;
     static float a3_oneemp = 1e10f;
     int condition;
     float tempa2, tempa3;
     float empinterp;

     if (whichemp == 0){
          condition = ( (a2_oldemp != a[2]) || (a3_oldemp != a[3]) || first );
     }
     if (whichemp == 1){
          condition = ( (a2_oneemp != a[2]) || (a3_oneemp != a[3]) || first );
     }

     if (condition){ 
          if (funny){
               fprintf(logfile, "a3 & a4 = %f %f \n", a2_oneemp, a3_oneemp);
               fprintf(logfile, "a(3) & a(4) = %f %f \n", a[2], a[3]);
               fprintf(logfile, "xy(1) & xy(2) = %d %d \n", xy[0], xy[1]);
          }
          first = 0; //false
          second = 1; //true
       
          if ( (*fitcall_ptr == 1) &&
               ((fabsf(a[2]) >= 4.0f) || (fabsf(a[3]) >= 4.0f)) ) {
               printf("improve proposes to move centroid significantly\n");
               printf("by (%f, %f), which would exceed memory ",a[2],a[3]);
               printf("allocations for empirical matrices\n");
               printf("artificially tempering fit to max move of 4.0\n");
               printf("strongly recommend limiting ABSLIM2 &3 in parameter ");
               printf("file and rerunning dophot\n");
               if (fabsf(a[2]) > 4.0f){
                    if (a[2] > 0.0f){
                         a[2] = 4.0f;
                    }
                    else{
                         a[2] = -4.0f;
                    }
               }
               if (fabsf(a[3]) > 4.0f){
                    if (a[3] > 0.0f){
                         a[3] = 4.0f;
                    }
                    else{
                         a[3] = -4.0f;
                    }
               }
          }

          p  = modff(-a[2], &flt_ip); 
          ip = (int)flt_ip;
          q  = modff(-a[3], &flt_iq); 
          iq = (int)flt_iq;
          if (funny){
               fprintf(logfile, "p, q, ip, & iq = %f %f %d %d \n", p, q, ip, iq);
          }
          if (p < 0.0f){
               p  += 1.0f;
               ip -= 1   ;
          }
          if (q < 0.0f){
               q  += 1.0f;
               iq -= 1   ;
          }
          /* fourpoint coefficients */
          c00 = (1.0f - p)*(1.0f - q);
          c01 = (1.0f - p)*(       q);
          c10 = (       p)*(1.0f - q);
          c11 = (       p)*(       q);
          /* handy defns for symmetry */
          q0 = q;
          p0 = p;
          q1 = 1.0f - q;
          p1 = 1.0f - p;
          /* twelveempoint coefficients */
          e10 = c00*q0*(q0 - 1.0f)/2.0f;
          e01 = c00*p0*(p0 - 1.0f)/2.0f;
          e23 = c11*q1*(q1 - 1.0f)/2.0f;
          e32 = c11*p1*(p1 - 1.0f)/2.0f;
          e20 = c10*q0*(q0 - 1.0f)/2.0f;
          e31 = c10*p1*(p1 - 1.0f)/2.0f;
          e13 = c01*q1*(q1 - 1.0f)/2.0f;
          e02 = c01*p0*(p0 - 1.0f)/2.0f;
          e11 = c00*(1.0f + p0*q0 - p0*p0 - q0*q0)
              + c11*(p1*q1)
              + c01*q1*(q1 - 2.0f*p0 + 1.0f)/2.0f
              + c10*p1*(p1 - 2.0f*q0 + 1.0f)/2.0f;
          e22 = c00*(p0*q0)
              + c11*(1.0f + p1*q1 - p1*p1 - q1*q1)
              + c01*p0*(p0 - 2.0f*q1 + 1.0f)/2.0f
              + c10*q0*(q0 - 2.0f*p1 + 1.0f)/2.0f;
          e12 = c00*q0*(q0 - 2.0f*p0 + 1.0f)/2.0f
              + c11*p1*(p1 - 2.0f*q1 + 1.0f)/2.0f
              + c01*(1.0f + p0*q1 - p0*p0 - q1*q1)
              + c10*p1*q0;
          e21 = c00*p0*(p0 - 2.0f*q0 + 1.0f)/2.0f
              + c11*q1*(q1 - 2.0f*p1 + 1.0f)/2.0f
              + c01*p0*q1
              + c10*(1.0f + p1*q0 - p1*p1 - q0*q0);

          tempa2 = a[2];
          tempa3 = a[3];
          if (whichemp == 0){
               a2_oldemp = tempa2;
               a3_oldemp = tempa3;
          }
          if (whichemp == 1){
               a2_oneemp = tempa2;
               a3_oneemp = tempa3;
          }

          r  = p - 0.5f;
          s  = q - 0.5f;
          ir = (int)r  ;
          is = (int)s  ;
          if (r < 0.0f){
               r  += 1.0f;
               ir -= 1   ;
          }
          if (s < 0.0f){
               s  += 1.0f;
               is -= 1   ;
          }
          d00 = (1.0f - r)*(1.0f - s);
          d01 = (1.0f - r)*(       s);
          d10 = (       r)*(1.0f - s);
          d11 = (       r)*(       s);
     }    

     ix1 = xy[0] + ip - 1 + (IHSIDE + 1);
     iy1 = xy[1] + iq - 1 + (IHSIDE + 1);

     ix0 = ix1 - 1;
     ix2 = ix1 + 1;
     ix3 = ix1 + 2;
     iy0 = iy1 - 1;
     iy2 = iy1 + 1;
     iy3 = iy1 + 2;
//     if( (ix0 <= 0) ||
//         (iy0 <= 0) ||
//         (ix3 >= (2*IHSIDE)) ||
//         (iy3 >= (2*IHSIDE)) ){ 
//          printf("x, y, from chisq = %d %d \n", xy[0], xy[1]);
//          printf("ip, iq,          = %d %d \n", ip   , iq);
//          printf("max, min x,y are %d %d %d %d \n", ix0, iy0, ix3, iy3);
//     }

     fa[1] = e11*(float)(mp[iy1][ix1])
           + e12*(float)(mp[iy2][ix1])
           + e21*(float)(mp[iy1][ix2])
           + e22*(float)(mp[iy2][ix2])
           + e01*(float)(mp[iy1][ix0])
           + e10*(float)(mp[iy0][ix1])
           + e02*(float)(mp[iy2][ix0])
           + e13*(float)(mp[iy3][ix1])
           + e20*(float)(mp[iy0][ix2])
           + e31*(float)(mp[iy1][ix3])
           + e32*(float)(mp[iy2][ix3])
           + e23*(float)(mp[iy3][ix2]);

     if (NEEDIT){
          iu0 = ix1 + ir;
          iv0 = iy1 + is;
          iu1 = iu0 + 1;
          iv1 = iv0 + 1;
          fa[2] = d00*(dxmp[iv0][iu0])
                + d01*(dxmp[iv1][iu0])
                + d10*(dxmp[iv0][iu1])
                + d11*(dxmp[iv1][iu1]);
          fa[3] = d00*(dymp[iv0][iu0])
                + d01*(dymp[iv1][iu0])
                + d10*(dymp[iv0][iu1])
                + d11*(dymp[iv1][iu1]);
          fa[2] = -a[1]*fa[2];
          fa[3] = -a[1]*fa[3];
          fa[0] = 1.0f;
     }
     
     empinterp = a[0] + a[1]*fa[1];

     if (funny && second){
          fprintf(logfile, "p0, q0, p1, q1 = %f %f %f %f \n", p0, q0, p1, q1); 
          fprintf(logfile, "r, s, ir, is = %f %f %d %d \n", r, s, ir, is); 
          fprintf(logfile, "the e coefficients: \n");
          fprintf(logfile, "%f %f %f %f \n", e10, e01, e23, e32); 
          fprintf(logfile, "%f %f %f %f \n", e31, e13, e02, e20); 
          fprintf(logfile, "%f %f %f %f \n", e11, e22, e12, e21); 
          fprintf(logfile, "the c coefficients: %f %f %f %f \n"
                            , c00, c01, c10, c11);
          fprintf(logfile, "ix0-3: %d %d %d %d \n", ix0, ix1, ix2, ix3); 
          fprintf(logfile, "iy0-3: %d %d %d %d \n", iy0, iy1, iy2, iy3); 
          second = 0; //false
     }

     /* reassigning common block vars */
     /* for oldemp */
     if (whichemp == 0){
          opqfake_.ip = ip;
          opqfake_.iq = iq;
          opqfake_.p  = p ;
          opqfake_.q  = q ;

          oefake_.e10 = e10;
          oefake_.e01 = e01;
          oefake_.e23 = e23;
          oefake_.e32 = e32;
          oefake_.e31 = e31;
          oefake_.e13 = e13;
          oefake_.e02 = e02;
          oefake_.e20 = e20;
          oefake_.e11 = e11;
          oefake_.e22 = e22;
          oefake_.e12 = e12;
          oefake_.e21 = e21;

          orsfake_.ir = ir;
          orsfake_.is = is;
          orsfake_.r  = r ;
          orsfake_.s  = s ;

          odfake_.d00 = d00;
          odfake_.d01 = d01;
          odfake_.d10 = d10;
          odfake_.d11 = d11;
     }
     /* for oneemp */
     if (whichemp == 1){
          pqfake_.ip = ip;
          pqfake_.iq = iq;
          pqfake_.p  = p ;
          pqfake_.q  = q ;

          efake_.e10 = e10;
          efake_.e01 = e01;
          efake_.e23 = e23;
          efake_.e32 = e32;
          efake_.e31 = e31;
          efake_.e13 = e13;
          efake_.e02 = e02;
          efake_.e20 = e20;
          efake_.e11 = e11;
          efake_.e22 = e22;
          efake_.e12 = e12;
          efake_.e21 = e21;

          rsfake_.ir = ir;
          rsfake_.is = is;
          rsfake_.r  = r ;
          rsfake_.s  = s ;

          dfake_.d00 = d00;
          dfake_.d01 = d01;
          dfake_.d10 = d10;
          dfake_.d11 = d11;
     }

     return empinterp;
}

double oldemp_(short int* xy, float* a, float* fa, int* m_ptr, int* fitcall_ptr)
{
     return empiricalinterp_(xy, a, fa, m_ptr, fitcall_ptr, 0);
}

double oneemp_(short int* xy, float* a, float* fa, int* m_ptr, int* fitcall_ptr)
{
     return empiricalinterp_(xy, a, fa, m_ptr, fitcall_ptr, 1);
}


