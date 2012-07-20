#ifndef PARUPD_H_
#define PARUPD_H_

/* dophot subroutine converted to c void function 02-05-2012 */
/*
C:  THIS SUBROUTINE IS PROFILE SPECIFIC AND NEEDS REWRITING FOR EACH
C:  FITTING FUNCTION.
*/

void parupd_(float* A, float* STARPAR, int* IX_ptr, int* IY_ptr);
void twoupd_(float* A, float* STARPAR, int* IX_ptr, int* IY_ptr);

#endif
