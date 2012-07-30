#include <stdlib.h>
#include "logh.h"
#include "tuneable.h"
#include "subraster_struct.h"
#include "crudestat_struct.h"
#include "unitize_struct.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "fillerup.h"

/* dophot subroutine converted to c void function 02-14-2012 */

void fillerup_(int** BIG, int** NOISE, int* IXIN_ptr, int* IYIN_ptr, int* JFAST_ptr, int* JSLOW_ptr)
{

     /* dereference pointers */
     int IXIN  = *IXIN_ptr ;
     int IYIN  = *IYIN_ptr ;
     int JFAST = *JFAST_ptr;
     int JSLOW = *JSLOW_ptr;

     /* renaming used commonblock variables in this function */
     short int* IRECT = tune2_.irect;
     float UFACTOR = unitize_.ufactor;
     int   NPT    = crudestat_.npt;
     float SUM0   = crudestat_.sum0;
     float SUM1   = crudestat_.sum1;
     float SUM2   = crudestat_.sum2;
     int   MAXVAL = crudestat_.maxval;
     float XM     = crudestat_.xm;
     float YM     = crudestat_.ym;
     float* Z  = subraster_.z ;
     float* YE = subraster_.ye;
     short int** XX = subraster_.xx;

     /* substance of subroutine begins here */
     short int IX, IY;

     float UFACTOR2;
     int J, JLO, JHI, I, ILO, IHI;
     int IBIG, INOISE;
     int IXI, IYI;

     UFACTOR2 = UFACTOR*UFACTOR;
     XM = 0.0f;
     YM = 0.0f;
     NPT = 0;
     SUM0 = 0.0f;
     SUM1 = 0.0f;
     SUM2 = 0.0f;

     JLO = max( IYIN - (int)(IRECT[1]/2) , 1    );
     JHI = min( IYIN + (int)(IRECT[1]/2) , JSLOW);
     ILO = max( IXIN - (int)(IRECT[0]/2) , 1    );
     IHI = min( IXIN + (int)(IRECT[0]/2) , JFAST);
     MAXVAL   = (int)BIG[JLO-1][ILO-1];
     for (J = JLO; J <= JHI; J++){
          for (I = ILO; I <= IHI; I++){
               IBIG   = (int)BIG[J-1][I-1];
               INOISE = NOISE[J-1][I-1];
               if (INOISE < MAGIC){
                    IX = (short int)(I - IXIN);
                    IY = (short int)(J - IYIN);
                    if (IBIG > MAXVAL){
                         IXI = (int)(IX);
                         IYI = (int)(IY);
                         if ( max( abs(IXI), abs(IYI) ) <= 1 ){
                              MAXVAL = IBIG;
                              XM     = (float)IX; //position of largest val
                              YM     = (float)IY; //position of largest val
                         }
                    }
                    NPT += 1;
                    XX[NPT-1][0] = IX;
                    XX[NPT-1][1] = IY;
                    Z[NPT-1]  = (float)(IBIG) / UFACTOR;
                    YE[NPT-1] = (float)(INOISE) / UFACTOR2;
                    SUM0 += 1.0f/YE[NPT-1];
                    SUM1 += (float)(IBIG) / YE[NPT-1];
               }
          }// end of I loop
     }// end of J loop

     if (SUM0 != 0.0f){
          SUM2 = SUM1/SUM0;
     }

     /* repointing changed common block variables and passed pointers */
     crudestat_.maxval = MAXVAL;
     crudestat_.xm   = XM;
     crudestat_.ym   = YM;
     crudestat_.npt  = NPT;
     crudestat_.sum0 = SUM0;
     crudestat_.sum1 = SUM1;
     crudestat_.sum2 = SUM2;

}

