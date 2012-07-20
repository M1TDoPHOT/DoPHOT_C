#ifndef ADDSTAR_H_
#define ADDSTAR_H_

/* dophot subroutine converted to c void fucntion 02-26-2011 */

// int output_img is a toggle of whether or not to output an image of the model
// subtracted or added to a file specified by img_file

void addstar_(double (*ONESTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, float* STARPAR, short int* ADDAREA, int* IADD_ptr, int model_img, char* model_file, int clean_img, char* clean_file);

#endif
