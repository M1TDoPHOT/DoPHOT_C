#ifndef ONEFIT_H_
#define ONEFIT_H_

double onefit_( double (*FUNCTN)(short int*, float*, float*, int*, int*), short int** X, float* Y, float* YE, int* N_ptr, float* A, float* FA, float* C_ptr, int* M_ptr, float* ACC, float* ALIM, int* IT_ptr, int which_model);

/* a wrapper function for chisq that will convert parameters and limits to
 and from log space if desired for fitting purposes */
// currently, update for each model switch, latter make toggle with flag[0]

#endif

