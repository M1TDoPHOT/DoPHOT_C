#ifndef WARMSTART_H_
#define WARMSTART_H_

/* dophot subroutine converted to c void functoin 03-30-2012 */

// WARNING not yet tested with a shadow file

//onestar is a pointer to the function name ONESTAR which returns double
void warmstart_( double (*ONESTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, char* filename );

#endif
