#ifndef OUTPUTS_H_
#define OUTPUTS_H_

/* dophot subroutine converted to c void function 02-06-2012 */

char* stdotpt_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPARR_ptr, float* APPARR, int* NAPPLE_ptr, float* PROBG_ptr, int which_model);

char* badotpt_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPMAX_ptr, float* APPARR, int* NAPPLE_ptr, float* PROBG_ptr);

char* sumout_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPARR_ptr, float* APPARR, int* NAPPLE_ptr, float* PROBG_ptr, int which_model);

char* shdout_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPARR_ptr, int which_model);

char* covarout_(float* SHADCOVAR);

#endif
