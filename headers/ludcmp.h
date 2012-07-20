#ifndef LUDCMP_H_
#define LUDCMP_H_

/* dophot subroutine converted to c void function 01-19-2012 */
/* subroutine used for matrix inversion, with lubksb.
   content lifted from numerical recipes in c/fortran        */ 

void ludcmp_(float **A, int N, int* INDX, float **lu);

#endif
