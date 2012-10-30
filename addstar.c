#include <stdlib.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "free_parking_struct.h"
#include "unitize_struct.h"
#include "drfake_struct.h"
#include "empmom_struct.h"
#include "changed_struct.h"
#include "funnypass_struct.h"
#include "galpass_struct.h"
#include "guess.h"
#include "addlims.h"
#include "empiricals.h" //contains oldemp
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "newfits.h"
#include "addstar.h"

/* dophot subroutine converted to c void fucntion 02-26-2011 */

// int output_img is a toggle of whether or not to output an image of the model
// subtracted or added to a file specified by img_file

void addstar_(double (*ONESTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int NFAST, int NSLOW, float* STARPAR, short int* ADDAREA, int IADD, int model_img, char* model_file, int clean_img, char* clean_file)
{

     /* rename used common block variables */
     int    NPAR = tune4_.npar;
     float  EFAC = tune22_.efac;
     float  GFAC = tune22_.gfac;
     float   FAC =  tune5_.fac;
     float  XPND =  tune5_.xpnd;
     float GXPND =  tune5_.gxpnd;
     float EPERDN = tune11_.eperdn;
     int  USEOLD = changed_.useold;
     int  NEEDIT = drfake_.needit;
     int BIGFOOT = galpass_.bigfoot;
     float E5 = empmom_.e5;
     float E6 = empmom_.e6;
     float E7 = empmom_.e7;
     float E2 = empmom_.e2;
     float E2OLD = empmom_.e2old;
     float UFACTOR = unitize_.ufactor;

     /* substance of subroutine begins here */
     short int* IX    = free_parking_.sifourarray_1;
     float* A    = free_parking_.npmaxarray_1;
     float* FA   = free_parking_.npmaxarray_2;
     float* B    = free_parking_.npmaxarray_3;

     int BADNEWS;
     int ZERO_dum = 0; //so can be passed by pointer to fxns
     int M; //so can be passed by pointer to fxns
     static int frst = 1; //true
     double guess_return;
     float SKY;
     int IXIN, IYIN;
     int KADD;
     float RFAC, TXPND;
     float BFACTOR, BSKY;
     float CFACTOR, CSKY;
     float E2C, TVAL, BTVAL;
     int JRECT0, JRECT1, JRECT2, JRECT3;
     int J, JHI, JLO; 
     int I, IHI, ILO; 
     int ITNOISE, ITEMP, IBIG, INOISE, IVAL, K;

     //for outputting added arrays to file if requested
     int z_x, z_y;
     int nx, ny;
     int** model_arr;
     int** clean_arr;
    
     if (frst){ 
          USEOLD = 0; //false
          frst = 0; //no longer first pass
     }
     NEEDIT = 0; //false

     /* A bit dangerous since a,b(5-7) aren't defined for empirical */
     guess_return = guess2_(A, STARPAR, &IXIN, &IYIN);
     guess_return = guess2_(B, STARPAR, &IXIN, &IYIN);
     SKY = (float)(guess_return);

     KADD = IADD/abs(IADD);
     /* We're going to save and recall limits;  seems devilishly simple;
     am I sure the parameters aren't corrupted?  IADD < 0 implies
     new subtraction, compute new limits.  IADD > 0 implies use
     old limits (and old function)! */
     if (abs(IADD) >= 2){
          RFAC = EFAC;
          if (IADD < 0){
               FA[0] = STARPAR[0]; 
               FA[1] = STARPAR[1]*10000.0f; 
               FA[2] = STARPAR[2]; 
               FA[3] = STARPAR[3]; 
               FA[4] = E5;
               FA[5] = E6;
               FA[6] = E7;
               addlims_(FA, ADDAREA);
          }
     }
     else{
          if (BIGFOOT){
               RFAC = GFAC;
               addlims_(STARPAR, ADDAREA);
          }
          else{
               RFAC = FAC;
               addlims_(STARPAR, ADDAREA);
          }
     } // end iadd if/else
               
     /* No changes here due to reordering of parameters.
        unused if empirical so not to worry */
     if (BIGFOOT){
          TXPND = GXPND;
     }
     else{
          TXPND = XPND;
     }

     B[4] = A[4]*TXPND*TXPND;
     B[5] = A[5]/(TXPND*TXPND);
     B[6] = A[6]*TXPND*TXPND;
     BFACTOR = UFACTOR*(float)KADD;
     BSKY    = (float)(KADD)*(-1.0f*SKY + 0.5f);
     CFACTOR = UFACTOR*RFAC;
     CSKY    = RFAC*(-1.0f*SKY) + 0.5f;
     JRECT0 = (int)ADDAREA[0];
     JRECT1 = (int)ADDAREA[1];
     JRECT2 = (int)ADDAREA[2];
     JRECT3 = (int)ADDAREA[3];

     ILO = max(JRECT0, 1); 
     IHI = min(JRECT1, NFAST); 
     JLO = max(JRECT2, 1); 
     JHI = min(JRECT3, NSLOW); 

     // for image output if requested
     nx = IHI - ILO + 1;
     ny = JHI - JLO + 1;
     if (model_img == 1){
          model_arr = malloc_int_2darr(ny, nx);
     }
     if (clean_img == 1){
          clean_arr = malloc_int_2darr(ny, nx);
     }

     for(J = JLO; J <= JHI; J++){
          z_y = J - JLO; // index counting from zero
          IX[1] = J - IYIN;
          for(I = ILO; I <= IHI; I++){
               z_x = I - ILO; // index counting from zero
               // in subtraction mode IBIG is the neighbor cleaned object
               IBIG   = (int)BIG[J-1][I-1];
               INOISE = NOISE[J-1][I-1];
               if (clean_img == 1){
                    clean_arr[z_y][z_x] = IBIG;
               }
               if (INOISE < MAGIC){
                    IX[0] = I - IXIN;
                    if (USEOLD){
                         E2C  = E2OLD;
                         TVAL = (float)(oldemp_(IX, A, FA, 
                                                &M, &ZERO_dum));
                    }
                    else{
                         E2C  = E2;
                         TVAL = (float)( (*ONESTAR)(IX, A, FA, 
                                                    &M, &ZERO_dum));
                    }
                    /* problem here: old e2 or new e2! */
                    if (abs(IADD) >= 2){
                         ITNOISE = (int)(
                                    (UFACTOR*TVAL*10000.0f - (SKY+0.5f))
                                   *(1.0f/EPERDN)
                                   *(STARPAR[1]/E2C)
                                   );
                         ITNOISE = max(ITNOISE, 0);
                    } 
                    else{ 
                         ITNOISE = 0;
                    }

                    // IVAL is what is added to img = -model
                    // ITEMP is image - model = residual
                    IVAL  = (int)(BFACTOR*TVAL + BSKY);  
                    ITEMP = IBIG + IVAL;
                    if (model_img == 1){
                         model_arr[z_y][z_x] = (-IVAL);
                    }
                    //from short int days  BADNEWS = ( (ITEMP >= 32767) || (ITEMP <= -32768) );
                    /* this segment probably not needed at all now with ints,
                       but is kept for historical reasons (laziness)
                       and for realistic check on data which should only go to
                       65535 realistically, although the int type may go to 
                       -2,147,483,648 to 2,147,483,647 at least in C.
                       BADNEWS will flag serious WTF errors */ 
                    BADNEWS = ( (ITEMP >= 66000) || (ITEMP <= -66000) );
                    if (BADNEWS){
                         fprintf(logfile,"I, J, & BADNEWS = %d %d %d ",
                                                   I, J, ITEMP);
                         fprintf(logfile,"(.....routine ADDSTAR)\n");
                         fprintf(logfile,"ix(1&2) = %d %d \n", IX[0], IX[1]);
                         fprintf(logfile,"a(n), n = 1,NPAR: ");
                         for(K = 0; K < NPAR; K++){
                              fprintf(logfile,"%f ", A[K]);
                         }
                         fprintf(logfile,"\n");
                         fprintf(logfile,"SERIOUS problem.. do not pass go ");
                         fprintf(logfile,"Do not collect 200\n");
                         fprintf(logfile,"kadd & bfactor = %d %f \n", 
                                                 KADD, BFACTOR);
                         fprintf(logfile,"tval & big(i,j) = %f %d \n", 
                                                 TVAL, IBIG);
                         fprintf(logfile,"ival & bsky = %d %f \n", 
                                                 IVAL, BSKY);
                         fprintf(logfile,"big(i,j) & noise(i,j) = %d %d \n", 
                                                 IBIG, INOISE);
                         // short int remnant IBIG   = -32768;
                         IBIG = OBLITVAL;
                         BIG[J-1][I-1] = IBIG;
                         INOISE = MAGIC;
                         NOISE[J-1][I-1] = INOISE;
                    }
                    else{
                         IBIG   = ITEMP;
                         BIG[J-1][I-1] = IBIG;

                         if (abs(IADD) >= 2){
                              IVAL = (int)(CFACTOR*TVAL + CSKY); 
                         } 
                         else{ 
                              BTVAL = (float)( (*ONESTAR)(IX, B, FA, 
                                                    &M, &ZERO_dum));
                              IVAL = (int)(CFACTOR*BTVAL + CSKY); 
                         } 
                         INOISE -= KADD*(IVAL*IVAL + ITNOISE);
                         NOISE[J-1][I-1] = INOISE;

                         if (INOISE <= 0){
                              fprintf(logfile,"I, J, & NEGNOISE = %d %d %d ",
                                       I, J, INOISE);
                              fprintf(logfile,"(.....routine ADDSTAR)\n");
                              fprintf(logfile,"SERIOUS problem.. do not pass go ");
                              fprintf(logfile,"Do not collect 200\n");
                              fprintf(logfile,"\n");
                         }
                    }// end BADNEWS if/else
               }// end NOISE < MAGIC if
          }//end I loop
     }//end J loop
     NEEDIT = 1; //true
                         
     /* repointing changed pointers and common block vars*/
     changed_.useold = USEOLD;
     drfake_.needit = NEEDIT;

     /* free locally allocated memory */
     if (model_img == 1){
          newfits_(nx, ny, model_arr, model_file, 0, " ");
          free_int_2darr(ny, model_arr);
     }
     if (clean_img == 1){
          newfits_(nx, ny, clean_arr, clean_file, 0, " ");
          free_int_2darr(ny, clean_arr);
     }

}
