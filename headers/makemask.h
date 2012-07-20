#ifndef MAKEMASK_H_
#define MAKEMASK_H_

/* dophot subroutine converted to c void functoin 02-17-2012 */

//onestar is a pointer to the function name ONESTAR which returns double
//onestar is usually pseud2d
//void makemask_( double (*ONESTAR)(short int* IX, float* A, float* FA, int* M_ptr, int* MMAX_ptr) )
void makemask_( double (*ONESTAR)(short int*, float*, float*, int*, int*) );

#endif
