#include <stdlib.h>
#include "logh.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "medfil.h"

/* dophot subroutine converted to c void function 01-23-11 */

//.... 32767 has often been written as 32770-3 to prevent I*2 overflow 
//           in computation

void medfil_(int* NFAST_ptr, int* NSLOW_ptr, int** LPICT, int** NEWPICT, int* Z1_ptr, int* Z2_ptr, int* XHW_ptr, int* YHW_ptr, int* nrmax_ptr, int* ncmax_ptr, int* prec_ptr)
{
     
     /* dereferencing the pointers */
         /* rows and columns are backwards because 
            they are backwards in the fortran 
            version... a mistake there? */
     int NFAST = *NFAST_ptr ; //NFAST = NX in image = NCOLS
     int NSLOW = *NSLOW_ptr ; //NSLOW = NY in image = NROWS
     int Z1    = *Z1_ptr    ;
     int Z2    = *Z2_ptr    ;
     int XHW   = *XHW_ptr   ;
     int YHW   = *YHW_ptr   ;
     int prec  = *prec_ptr  ;

     /* substance of subroutine begins here */
     short SC_1 = 0; //flags for 'short circuits'
     short SC_2 = 0; //used in place of gotos
     short SC_3 = 0;
     int i, j, krow, kcol;
     int m, count ;
     int k, kc, klim     ;
     int min4 = 0 ; 
     short int* BOX = malloc_si_1darr(65536);

     Z1 = max(Z1,-32767) ;
     Z2 = min(Z2, 32767) ;

     for (i = YHW; i < (NSLOW-YHW); i++){

          if ((i != YHW) && ( (i % 50) == 0)){
               fprintf(logfile, "+    %7d   ->   %7d %7d %7d %7d %7d \n", 
                                   prec*min4, i, 
                                   NEWPICT[i - 1][(  NFAST/4) - 1],
                                   NEWPICT[i - 1][(  NFAST/2) - 1],
                                   NEWPICT[i - 1][(3*NFAST/4) - 1],
                                   min4) ;
          }

          /* Load fresh box : */
          /* box will store number of counts for each DN
             between max and min allowable */
          for (m = 0; m < 65535; m++){
               BOX[m] = 0 ;
          }

          min4 = Z2/prec ;
          count = 0 ;
          for (krow = (i - YHW); krow < (i + YHW + 1); krow++){
               for (kcol = 0; kcol < (2*XHW + 1); kcol++){
                    if ((LPICT[krow][kcol] >= Z1) &&
                        (LPICT[krow][kcol] <= Z2)  ){
                         BOX[(LPICT[krow][kcol]/prec) + 32767] += 1 ;
                         if (LPICT[krow][kcol]/prec <= min4){
                              min4 = LPICT[krow][kcol]/prec ;
                         }
                         count += 1 ;
                    }
               }
          }

          /* if box is empty, short-circuit */
          /* Set LOW so that objects are triggered 
                  on 1st pass but rejected on 
                  AVERAGE sky test.... setting high 
                  may have other problems   */
          SC_1 = 0 ;
          if (count == 0){
               NEWPICT[i][XHW] = 0 ;
               SC_1 = 1   ;
          }

          kc = 0; 
          k = (min4 + 32767) ;
          klim = count/2 + 1 ;
          while ( (SC_1 == 0) &&  (k < (Z2/prec + 32767 + 1)) ){
               kc += BOX[k] ;
               if (kc >= klim){
                    NEWPICT[i][XHW] = (k - 32767)*prec ;
                    SC_1 = 1 ;
               }
               k++ ;
          }

          for (j = (XHW + 1); j < (NFAST - XHW); j++){
               /* strip old left column and Add new right column */
               for (k = (i - YHW); k < (i + YHW + 1); k++){
                    if ( (LPICT[k][j - (XHW + 1)] >= Z1) && 
                         (LPICT[k][j - (XHW + 1)] <= Z2) ){
                         BOX[LPICT[k][j - (XHW + 1)]/prec + 32767] -= 1 ;
                         count -= 1 ;
                    }
                    if ( (LPICT[k][j + XHW] >= Z1) && 
                         (LPICT[k][j + XHW] <= Z2) ){
                         BOX[LPICT[k][j + XHW]/prec + 32767] += 1 ;
                         if (LPICT[k][j + XHW]/prec <= min4){
                              min4 = LPICT[k][j + XHW]/prec ;
                         }
                         count += 1 ;
                    }
               }

               SC_2 = 0 ;
               /* if box is empty, short-circuit */
               if (count == 0){
                    NEWPICT[i][j] = 0 ;
                    SC_2 = 1   ;
               }

               kc = 0;
               klim = count/2 + 1 ;
               k = (min4 + 32767) ;
               while ( (SC_2 == 0) && (k < (Z2/prec + 32767 + 1)) ){
                    kc += BOX[k] ;
                    if (kc >= klim){
                         NEWPICT[i][j] = (k - 32767)*prec ; 
                         SC_2 = 1 ;
                    }
                    k++ ; 
               }

          }// end of j loop
     }// end of i loop

     /* fill in the margins */
     for (i = 0; i < NSLOW; i++){
          for (j = 0; j < NFAST; j++){
               SC_3 = 0 ;
               if( (i >= YHW) && ( i <= NSLOW - YHW - 1) &&
                   (j >= XHW) && ( j <= NFAST - XHW - 1) ){
                    SC_3 = 1 ;
               }

               if ( (SC_3 == 0) && (i <= (YHW - 1)) ){
                    if ( (SC_3 == 0) && (j <= (XHW - 1)) ){
                         NEWPICT[i][j] = NEWPICT[YHW][XHW];
                         SC_3 = 1 ;
                    }
                    if ( (SC_3 == 0) && (j >= (NFAST - XHW - 1)) ){
                         NEWPICT[i][j] = NEWPICT[YHW][NFAST - XHW - 1];
                         SC_3 = 1 ;
                    }
                    if ( (SC_3 == 0) ){
                         NEWPICT[i][j] = NEWPICT[YHW][j] ; 
                         SC_3 = 1 ;
                    }
               }

               if ( (SC_3 == 0) && (i >= (NSLOW - YHW - 1)) ){
                    if ( (SC_3 == 0) && (j <= (XHW - 1)) ){
                         NEWPICT[i][j] = NEWPICT[NSLOW - YHW - 1][XHW];
                         SC_3 = 1 ;
                    }
                    if ( (SC_3 == 0) && (j >= (NFAST - XHW - 1)) ){
                         NEWPICT[i][j] = 
                              NEWPICT[NSLOW - YHW - 1][NFAST - XHW - 1];
                         SC_3 = 1 ;
                    }
                    if ( (SC_3 == 0) ){
                         NEWPICT[i][j] = NEWPICT[NSLOW - YHW - 1][j] ; 
                         SC_3 = 1 ;
                    }
               }

               if ( (SC_3 == 0) && (j <= (XHW - 1)) ){
                    NEWPICT[i][j] = NEWPICT[i][XHW] ; 
                    SC_3 = 1 ;
               }
               if ( (SC_3 == 0) && (j >= (NFAST - XHW - 1)) ){
                    NEWPICT[i][j] = NEWPICT[i][NFAST - XHW - 1] ; 
                    SC_3 = 1 ;
               }

          }// end j loop, reset SC_3
     }// end i loop
                            
     /* repointing the changed pointers */
     *Z1_ptr = Z1 ;
     *Z2_ptr = Z2 ;
     /* freeing newly allocated memory */
     free(BOX) ;

     return;
}
