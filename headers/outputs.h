#ifndef OUTPUTS_H_
#define OUTPUTS_H_

/* dophot subroutine converted to c void function 02-06-2012 */

char* stdotpt_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPARR_ptr, float* APPARR, int* NAPPLE_ptr, float* PROBG_ptr);

char* badotpt_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPMAX_ptr, float* APPARR, int* NAPPLE_ptr, float* PROBG_ptr);

char* sumout_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPARR_ptr, float* APPARR, int* NAPPLE_ptr, float* PROBG_ptr);

char* shdout_(int* I_ptr, int* ITYPE_ptr, float* STPARR, int* NPARR_ptr);

#endif
