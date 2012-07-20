#include <stdio.h>
#include "cast_arr.h"
/* dophot subroutine converted to c void function 01-18-2012 */

void lubksb_(float** A, int N, int* INDX, float* B, float* realB)
{

     /* substance of the suboutine begins here */
     int II,LL ;
     int I,J   ;
     float SUM    ;
     II = -1      ;    
 
     for (I = 0; I < N; I++){
          LL    = INDX[I] - 1;
          SUM   = B[LL]   ;
          B[LL] = B[I]    ;
          if (II != -1){
               for (J = II; J < I; J++){
                    SUM = SUM - A[J][I]*B[J] ;
               }
          } 
          else{
               if (SUM != 0.0f){
                    II = I ;    
               }
          }
          B[I]  = SUM     ;
     }

     for (I = (N-1); I > -1; I--){
          SUM   = B[I]    ;
          if (I < (N-1)){
               for (J = I+1; J < N; J++){
                    SUM = SUM - A[J][I]*B[J] ;
               }
          }
          B[I] = SUM / A[I][I] ;
     }
          
}

