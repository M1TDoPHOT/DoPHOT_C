#include <stdlib.h>
#include <math.h>
#include "tuneable.h"
#include "logh.h"
#include "cast_arr.h"
#include "ludcmp.h"

/* dophot subroutine converted to c void function 01-19-2012 */
/* subroutine used for matrix inversion, with lubksb.
   content lifted from numerical recipes in c/fortran        */ 


void ludcmp_(float **A, int N, int* INDX, float** lu)
{

     /* substance of subroutine begins here */
     int   I,J,K,IMAX    ;
     float AAMAX,SUM,DUM ;
     float* VV = malloc_float_1darr((N+1)) ;
     if (VV == NULL){
          fprintf(stderr, "out of memory\n") ;
     } 
     IMAX = -1   ;

     for (I = 0; I < N; I++){
          AAMAX = 0.0f ;
          for (J = 0; J < N; J++){
               if (fabsf(A[J][I]) > AAMAX){
                    AAMAX = fabsf(A[J][I]) ; 
               }
          }
          if (AAMAX == 0.0){
               fprintf(logfile, "Singular Matrix\n") ;
          }
          VV[I] = 1.0f/AAMAX ;
     }

     for (J = 0; J < N; J++){
          if (J > 0){
               for (I = 0; I < J; I++){
                    SUM = A[J][I] ;
                    if (I > 0){
                         for (K = 0; K < I; K++){
                              SUM = SUM - A[K][I]*A[J][K] ;
                         }
                         A[J][I] = SUM ;
                    }
               }
          }
          
          AAMAX = 0.0 ;
          for (I = J; I < N; I++){
               SUM = A[J][I] ;
               if (J > 0){
                    for (K = 0; K < J; K++){
                         SUM = SUM - A[K][I]*A[J][K] ;
                    }
                    A[J][I] = SUM ;
               } 
               DUM = VV[I]*fabsf(SUM) ;
               if (DUM >= AAMAX){
                    IMAX  = I   ;
                    AAMAX = DUM ;
               }
          }
//          printf("IMAX, AAMAX = %d, %f \n", IMAX, AAMAX) ;

          if (J != IMAX){
               for (K = 0; K < N; K++){
                    DUM        = A[K][IMAX] ;
                    A[K][IMAX] = A[K][J]    ;
                    A[K][J]    = DUM        ;
               }
               VV[IMAX] = VV[J] ;
          }

          // WARNING this line will only work so long as INDX is
          // only holding index positions in 1 conuting scheme
          INDX[J] = IMAX + 1 ;

          if (J != (N-1)){
               if (A[J][J] == 0.0f){
                    A[J][J] = TINY ;
               }
               DUM = 1.0f/A[J][J] ;
               for (I = (J+1); I < N; I++){
                    A[J][I] = A[J][I]*DUM ;
               }
          }
     }

     if (A[N-1][N-1] == 0.0f){
          A[N-1][N-1] = TINY ;
     }

     free(VV)   ;

}
