#include <stdlib.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "free_parking_struct.h"
#include "empmom_struct.h"
#include "guess.h"
#include "addlims.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "newfits.h"
#include "stampout.h"

/* for outputting postage stamps of the model size subtracted of each object */

// int output_img is a toggle of whether or not to output an image of the model
// subtracted or added to a file specified by img_file

void stampout_(int** BIG, int* NFAST_ptr, int* NSLOW_ptr, float* STARPAR, short int* ADDAREA, int ITYPE, char* out_file)
{
     /* dereference pointers */
     int NFAST = *NFAST_ptr;
     int NSLOW = *NSLOW_ptr;
     /* rename used common block variables */
     float  GFAC = tune22_.gfac;
     float   FAC =  tune5_.fac;
     float E5 = empmom_.e5;
     float E6 = empmom_.e6;
     float E7 = empmom_.e7;

     /* substance of subroutine begins here */
     float* FA   = free_parking_.npmaxarray_1;

     float RFAC;
     int JRECT0, JRECT1, JRECT2, JRECT3;
     int J, JHI, JLO; 
     int I, IHI, ILO; 
     int z_x, z_y;
     int nx, ny;
     int** out_arr;
    
     if (abs(ITYPE) == 10){ //empirical star
          FA[0] = STARPAR[0];
          FA[1] = STARPAR[1]*10000.0f;
          FA[2] = STARPAR[2];
          FA[3] = STARPAR[3];
          FA[4] = E5;
          FA[5] = E6;
          FA[6] = E7;
          addlims_(FA, ADDAREA);
     }
     else{ 
          if (ITYPE == 2) { //galaxy
               RFAC = GFAC;
               addlims_(STARPAR, ADDAREA);
          }
          else{ //not galaxy
               RFAC = FAC;
               addlims_(STARPAR, ADDAREA);
          }
     } // end iadd if/else
               
     JRECT0 = (int)ADDAREA[0];
     JRECT1 = (int)ADDAREA[1];
     JRECT2 = (int)ADDAREA[2];
     JRECT3 = (int)ADDAREA[3];

     ILO = max(JRECT0, 1); 
     IHI = min(JRECT1, NFAST); 
     JLO = max(JRECT2, 1); 
     JHI = min(JRECT3, NSLOW); 

     nx = IHI - ILO + 1;
     ny = JHI - JLO + 1;
     out_arr = malloc_int_2darr(ny, nx);

     for(J = JLO; J <= JHI; J++){
          z_y = J - JLO; // index counting from zero
          for(I = ILO; I <= IHI; I++){
               z_x = I - ILO; // index counting from zero
               out_arr[z_y][z_x] = BIG[J-1][I-1];
          }//end I loop
     }//end J loop
                         
     newfits_(nx, ny, out_arr, out_file, 0, " ");

     /* free locally allocated memory */
     free_int_2darr(ny, out_arr);

}
