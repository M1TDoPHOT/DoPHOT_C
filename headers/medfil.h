#ifndef MEDFIL_H_
#define MEDFIL_H_

/* dophot subroutine converted to c void function 01-23-11 */

//.... 32767 has often been written as 32770-3 to prevent I*2 overflow 
//           in computation

void medfil_(int* NSLOW_ptr, int* NFAST_ptr, int** LPICT, int** NEWPICT, int* Z1_ptr, int* Z2_ptr, int* XHW_ptr, int* YHW_ptr, int* nrmax_ptr, int* ncmax_ptr, int* prec_ptr);

#endif
