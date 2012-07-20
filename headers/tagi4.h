#ifndef TAGI4_H_
#define TAGI4_H_

/* dophot subroutine converted to c void function 01-20-2012 */
/* SORTS ARRAY A INTO INCREASING ORDER, FROM A(II) TO A(JJ)
   ARRAYS IU(K) AND IL(K) PERMIT SORTING UP TO 2**(K+1)-1 ELEMENTS
   ARRAY TAG IS PERMUTED IN THE SAME MANNER AS ARRAY A  */

int compare (const void * qwerty, const void * poiuyt);

void tagi4_(int* A, int* II_ptr, int* JJ_ptr, short int* TAG);

#endif
