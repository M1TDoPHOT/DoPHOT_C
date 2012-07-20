#ifndef SHAPE_H_
#define SHAPE_H_

/* dophot subroutine converted to c void function 3-13-2012 */
//onestar is a pointer to the function name ONESTAR which returns double
//onestar is usually pseud2d
//twostar is usually pseud4d which takes the same parameters
void shape_( double (*ONESTAR)(short int*, float*, float*, int*, int*), double (*TWOSTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr );

#endif
