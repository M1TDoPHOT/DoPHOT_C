#ifndef COSMIC_H_
#define COSMIC_H_

/* dophot logical function converted to c int function 02-21-2012 */
//onestar is a pointer to the function name ONESTAR which returns double
//onestar is usually pseud2d
int cosmic_( double (*ONESTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, float* STARPAR );

#endif
     

