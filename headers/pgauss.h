#ifndef PGAUSS_H_
#define PGAUSS_H_

/* dophot function with 'entry' converted to 2 c double functions
   02-17-2012 */

double pgauss2d_(short int* IX, float* A, float* FA, int* M_ptr, int* fitcall_ptr);

double pgauss4d_(short int* IX, float* A, float* FA, int* M_ptr, int* fitcall_ptr);

#endif
