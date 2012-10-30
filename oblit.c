#include <stdlib.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "free_parking_struct.h"
#include "mini_mathlib.h"
#include "cast_arr.h"
#include "guess.h"
#include "ellipint.h"
#include "oblit.h"

/* dophot logical function converted to c int function 02-20-2012 */

//onestar is a pointer to the function name ONESTAR which returns double
//onestar is usually pseud2d
int oblit_( double (*ONESTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, float* STARPAR )
{

     /* dereference pointers */
     int NFAST = *NFAST_ptr;
     int NSLOW = *NSLOW_ptr;
     /* rename used common block variables */
     int lverb = tune14_.lverb;

     /* substance of subroutine begins here */
     float* A  = free_parking_.npmaxarray_1;
     int OBLIT;
     double guess2_return;
     int IX, IY;
     float WX, WY;
     int IXHI, IXLO, IYHI, IYLO;
     int JX, JY;

     guess2_return = guess2_(A, STARPAR, &IX, &IY);
     WX = STARPAR[4];
     WY = STARPAR[6];
     IXHI = min( (int)(IX + WX/2.0f + 0.5f) , NFAST);
     IXLO = max( (int)(IX - WX/2.0f + 0.5f) , 1    );
     IYHI = min( (int)(IY + WY/2.0f + 0.5f) , NSLOW);
     IYLO = max( (int)(IY - WY/2.0f + 0.5f) , 1    );
    
     if (lverb > 30){
          fprintf(logfile, "Obliterating following region :\n");
          fprintf(logfile, "IXLO, IXHI, IYLO, IYHI = %d, %d, %d, %d\n",
                            IXLO, IXHI, IYLO, IYHI);
     }

     for (JY = IYLO-1; JY < IYHI; JY++){
          for (JX = IXLO-1; JX < IXHI; JX++){
//               BIG[JY][JX] = -32768; 
               BIG[JY][JX] = OBLITVAL; 
               NOISE[JY][JX] = MAGIC; 
          }
     }

     if (lverb > 40){
          fprintf(logfile, "OBLITERATION: IX & IY = %d %d\n", IX, IY);
          fprintf(logfile, "OBLITERATION: WX & WY = %f %f\n", WX, WY);
     }
 
     OBLIT = 1; //true

     return OBLIT;
}


int oblit_ellipse_( double (*ONESTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, float* STARPAR )
{

     /* dereference pointers */
     int NFAST = *NFAST_ptr;
     int NSLOW = *NSLOW_ptr;
     /* rename used common block variables */
     int lverb = tune14_.lverb;

     /* substance of subroutine begins here */
     float* A  = free_parking_.npmaxarray_1;
     int OBLIT;
     double guess2_return;
     int IX, IY;
     int JX, JY;
     int IXHI, IXLO, IYHI, IYLO;
     float a, b, tilt, area;
     float sigmax, sigmay, sigmaxy;
     float dist;

     guess2_return = guess2_(A, STARPAR, &IX, &IY);
     a    = STARPAR[4];
     b    = STARPAR[6];
     tilt = STARPAR[5];
    
     ellipint_(&a, &b, &tilt, &area, &sigmax, &sigmaxy, &sigmay);
     sigmax = 1.0f/sigmax;
     sigmay = 1.0f/sigmay;

     // allowing a search region of 1.5x the x,y dimensions of the 
     // semi-major and minor axes
     IXHI = min( (int)(IX + (2.5f*a + 0.5f)) , NFAST);
     IXLO = max( (int)(IX - (2.5f*a + 0.5f)) , 1    );
     IYHI = min( (int)(IY + (2.5f*b + 0.5f)) , NSLOW);
     IYLO = max( (int)(IY - (2.5f*b + 0.5f)) , 1    );

     if (lverb > 30){
          fprintf(logfile, "Obliterating an ellipse within the following region :\n");
          fprintf(logfile, "IXLO, IXHI, IYLO, IYHI = %d, %d, %d, %d\n",
                            IXLO, IXHI, IYLO, IYHI);
     }

     for (JY = IYLO-1; JY < IYHI; JY++){
          for (JX = IXLO-1; JX < IXHI; JX++){
               dist = 0.5f*(  (JX - IX)*(JX - IX)*sigmax 
                            + 2.0f*sigmaxy*(JX - IX)*(JY - IY)
                            + (JY - IY)*(JY - IY)*sigmay);
               if (dist <= 1.0f){
                    BIG[JY][JX] = OBLITVAL;
                    NOISE[JY][JX] = MAGIC;
               }
          }
     }


     if (lverb > 40){
          fprintf(logfile, "OBLITERATION: IX & IY = %d %d\n", IX, IY);
          fprintf(logfile, "OBLITERATION: WX, WY, TILT = %f %f %f\n", a, b, tilt);
     }
 
     OBLIT = 1; //true

     return OBLIT;

}
