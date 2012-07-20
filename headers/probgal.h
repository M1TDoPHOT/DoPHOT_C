#ifndef PROBGAL_H_
#define PROBGAL_H_

/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:I THINK I WANT TO CALCULATE THE DERIVATIVE OF THE BEST FITTING CHISQUARED
C:WITH RESPECT TO A VARIATION IN A(2),A(5),A(6) AND A(7) WHICH KEEPS FLUX
C:CONSTANT.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

/* dophot function converted to c double function 02-20-2012 */
//onestar is a pointer to the function name ONESTAR which returns double
//onestar is usually pseud2d

double probgal_(double (*ONESTAR)(short int*, float*, float*, int*, int*), short int** XX, float* YY, float* ZZ, int* NPT_ptr, float* A, float* FA);

#endif
