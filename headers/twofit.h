#ifndef TWOFIT_H_
#define TWOFIT_H_

/* dophot function converted to c void functoin 02-22-2012 */

//twostar is a pointer to the function name TWOSTAR which returns double
//twostar is usually pseud4d
double twofit_( double (*TWOSTAR)(short int*, float*, float*, int*, int*), float* STARPAR );

#endif
