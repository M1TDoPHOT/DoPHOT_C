#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "starlist_struct.h"
#include "aperlist_struct.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "elarea.h"
#include "impaper.h"

/* dophot subroutine converted to c void function 01-29-2012 */
/* 
c
c  This is an independent way of calculating the difference between the
c  aperture and fit mags.  Large values for abs(APPLE[i][2]) for STAR 
c  objects will indicate problems.  This value equals the difference 
c  between the sky and fit mags.
c
*/

void impaper_(int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, int* K_ptr)
{

     /* dereference pointers */
     int NFAST = *NFAST_ptr; 
     int NSLOW = *NSLOW_ptr; 
     int K     = *K_ptr    ;
     /* shortening names of common block variables for easier use */
     float* ARECT = tune2_.arect;
     int    lverb = tune14_.lverb;
     float** STARPAR = starlist_.starpar;
     float** APPLE = aperlist_.apple;
     
     /* substance of subroutine begins here */
     double SUMAN, SUMBN, VARA, VARB;
     double STARLUM, RATLUM;
     float SUMA0, SUMA1 ;
     float SUMB0, SUMB1 ;
     int I, IHI, ILO, J, JHI;
     float TXA, TXB, TYA, TYB;
     float XA, XB, YA, YB, XYA, XYB;
     float ABJYC, ABIXC, DXY;
     int BADPIX = 0 ; //0 for false, 1 for true

     for (I = 0; I < NAPPLE; I++){
          APPLE[K-1][I] = 0;
     } 
     /* changed indices */
     float XC = STARPAR[K-1][2]; 
     float YC = STARPAR[K-1][3]; 

     SUMAN = 0.0;
     SUMBN = 0.0;
     SUMA0 = 0.0f;
     SUMA1 = 0.0f;
     SUMB0 = 0.0f;
     SUMB1 = 0.0f;

     JHI = min( (int)(YC + ARECT[1] + 0.5f), NSLOW );
     J   = max( (int)(YC - ARECT[1] + 0.5f), 1 ); 
     IHI = min( (int)(XC + ARECT[0] + 0.5f), NFAST );
     ILO = max( (int)(XC - ARECT[0] + 0.5f), 1 ); 

     TXA = ARECT[0]*0.5f + 0.5f ;
     TXB = ARECT[0]      + 0.5f ;
     TYA = ARECT[1]*0.5f + 0.5f ;
     TYB = ARECT[1]      + 0.5f ;

     while (J <= JHI){
          ABJYC = fabsf(J - YC);
          YA    = min1(1.0f, max1(0.0f, (TYA - ABJYC)) );
          YB    = min1(1.0f, max1(0.0f, (TYB - ABJYC)) );
          I     = ILO;
          while (I <= IHI){
               ABIXC = fabsf(I - XC);
               XA    = min1(1.0f, max1(0.0f, (TXA - ABIXC)) );
               XB    = min1(1.0f, max1(0.0f, (TXB - ABIXC)) );
               XYA   = XA*YA;
               XYB   = XB*YB;

               BADPIX = 0; // false unless true
               if ((NOISE[J-1][I-1] >= MAGIC) || 
                   (NOISE[J-1][I-1] <  0)      ){
                    BADPIX = 1; //true
               }

               if (XYA > 0){
                    if (BADPIX == 1){
                         I = IHI;
                         J = JHI;
                    }
                    else{
                         SUMA0 += XYA;
                         SUMA1 += XYA*BIG[J-1][I-1];
                         SUMAN += XYA*XYA*NOISE[J-1][I-1];
                    }
               }
               if (XYA < 1){
                    if (BADPIX == 0){
                         DXY    = (XYB - XYA)/NOISE[J-1][I-1];
                         SUMB0 += DXY;    
                         SUMB1 += DXY*BIG[J-1][I-1];    
                         SUMBN += DXY*DXY*NOISE[J-1][I-1];
                    }
               }
               I++ ;
          } // end of I while
          J++ ;
     }  // end of J while

     if (BADPIX != 1){
          if ( (fabsf( SUMA0/(ARECT[0]*ARECT[1]) - 1.0f)) > 0.001f){
               if (lverb > 20){
                    fprintf(logfile,"AREAS NOT EQUAL! STAR #, SUMB1 & SUMB0 = %5d %f %f\n",
                                    K, SUMB1, SUMB0);
               }
          } 
          else{
               APPLE[K-1][0] = SUMA1 - SUMB1*(SUMA0/SUMB0);
               APPLE[K-1][1] = SUMB1/SUMB0;
               VARA = SUMAN;
               VARB = (SUMA0/SUMB0)*(SUMA0/SUMB0)*SUMBN;
               APPLE[K-1][4] = 1.086f*(float)(sqrt(VARA + VARB))/APPLE[K-1][0];
               if (STARPAR[K-1][4] > 0.0f){
                    STARLUM  = elarea_(STARPAR[K-1]+4,STARPAR[K-1]+5,STARPAR[K-1]+6);
                    STARLUM  = 6.283185*STARLUM*STARPAR[K-1][1];
                    RATLUM   = STARLUM/APPLE[K-1][0];
                    if (RATLUM > 0.0){
                         APPLE[K-1][2] = 2.5f*log10f((float)RATLUM);
                    }
               }
          }
     } //end of if BADPIX ==0

     
}     
