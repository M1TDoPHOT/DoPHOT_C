#include <stdio.h>
#include <stdlib.h>
#include "cast_arr.h"
#include "tagi4.h"

/* dophot subroutine converted to c void function 01-20-2012 */
/* SORTS ARRAY A INTO INCREASING ORDER, FROM A(II) TO A(JJ)
   ARRAYS IU(K) AND IL(K) PERMIT SORTING UP TO 2**(K+1)-1 ELEMENTS
   ARRAY TAG IS PERMUTED IN THE SAME MANNER AS ARRAY A  */

int compare (const void * qwerty, const void * poiuyt);

void tagi4_(int* A, int* II_ptr, int* JJ_ptr, short int* TAG)
{

     /* dereferencing the non 1darray pointers */
     int II = *II_ptr ; 
     int JJ = *JJ_ptr ; 

//     printf("II, JJ = %d, %d, \n", II, JJ);
     /* ignoring the original dophot code and instead using 
        the c std library routine qsort */
     int size_of_a = (JJ-II+1);
     int* unsorted_A = malloc_int_1darr(size_of_a);
     if (unsorted_A == NULL){
          fprintf(stderr, "out of memory\n"); 
     }
     int* sorted_A = malloc_int_1darr(size_of_a);
     if (sorted_A == NULL){
          fprintf(stderr, "out of memory\n"); 
     }
     short int* unsorted_TAG = malloc_si_1darr(size_of_a);
     if (unsorted_TAG == NULL){
          fprintf(stderr, "out of memory\n"); 
     }
     short int* sorted_TAG = malloc_si_1darr(size_of_a);
     if (sorted_TAG == NULL){
          fprintf(stderr, "out of memory\n"); 
     }
     
     int i,j ;
     for (i = 0; i < size_of_a; i++){
          unsorted_A[i]   =   A[i + II - 1] ;
          sorted_A[i]     =   A[i + II - 1] ;
          unsorted_TAG[i] = TAG[i + II - 1] ;
     } 

     /* sort A */
     qsort(sorted_A, size_of_a, sizeof(int), compare) ;
     /* use old and new locationd of vals in A to sort TAG */
     /* by permuting through the new A list and matching 
        to values in the old, unsorted A list. */
     int A_val ;
     for (i = 0; i < size_of_a; i++){
          A_val = sorted_A[i] ; 
          for (j = 0; j < size_of_a; j++){
               if (A_val == unsorted_A[j]){
                    sorted_TAG[i] = unsorted_TAG[j] ;
                    A_val  = -9999         ;
               }
          }
          if (A_val != -9999){
               fprintf(stderr, "problem matching sort\n") ;
          }
     }

     /* repointing the passed pointers */
     for (i = 0; i < size_of_a; i++){
          TAG[i + II - 1] = sorted_TAG[i] ;
          A[i + II - 1]   = sorted_A[i]   ;
     }

     free(unsorted_A);
     free(sorted_A);
     free(unsorted_TAG);
     free(sorted_TAG);
     return ;   
}

int compare (const void * qwerty, const void * poiuyt)
{
  return ( *(int*)qwerty - *(int*)poiuyt );
}
