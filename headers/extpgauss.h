#ifndef EXTPGAUSS_H_
#define EXTPGAUSS_H_

/* dophot extended pseudogaussian function with b4 and b6 as free parameters 
   07-19-2012 */

double extpgauss2d_(short int* IX, float* A, float* FA, int* M_ptr, int* fitcall_ptr);

double extpgauss4d_(short int* IX, float* A, float* FA, int* M_ptr, int* fitcall_ptr);

#endif
