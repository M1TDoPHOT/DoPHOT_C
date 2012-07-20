#include "tuneable.h"
#include "offpic.h"

/* dophot function converted to c int function 02-02-2012 */
/* COMPUTES DISTANCES BEYOND EDGE OF PICTURE */

int offpic_(float* A, int* IX_ptr, int* IY_ptr, int* NFAST_ptr, int* NSLOW_ptr, float* DX_ptr, float* DY_ptr)
{

     /* dereferencing pointers */
     int IX    = *IX_ptr;
     int IY    = *IY_ptr;
     int NFAST = *NFAST_ptr;
     int NSLOW = *NSLOW_ptr;
     float DX  = *DX_ptr;
     float DY  = *DY_ptr;

     /* renaming commonblock variables */
     int fixpos = tune21_.fixpos;

     /* substance of function begins here */
     int NOGOOD;
     float X, Y;

     X = IX + A[2];
     Y = IY + A[3];

     if (X < 0.0f){
          DX = -X;
     }
     else{
          if (X > (float)NFAST){
               DX = X - NFAST;
          }
          else{
               DX = 0.0f;
          }
     }

     if (Y < 0.0f){
          DY = -Y;
     }
     else{
          if (Y > (float)NSLOW){
               DY = Y - NSLOW;
          }
          else{
               DY = 0.0f;
          }
     }

     NOGOOD = ( (DX != 0.0f) || (DY != 0.0f) );
     if (!(fixpos)){
          NOGOOD = ( NOGOOD || (A[1] < 0.0f) );
     }

     /* rereference changed pointers */
     *DX_ptr = DX;
     *DY_ptr = DY;

     return NOGOOD;
}
	
