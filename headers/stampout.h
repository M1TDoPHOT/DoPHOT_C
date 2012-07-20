#ifndef STAMPOUT_H_
#define STAMPOUT_H_

/* hack of addstar for outputting postage stamps of the model size subtracted of each object */

// int output_img is a toggle of whether or not to output an image of the model
// subtracted or added to a file specified by img_file
// IADD is >=2 for empirical stars, 1 else

void stampout_(int** BIG, int* NFAST_ptr, int* NSLOW_ptr, float* STARPAR, short int* ADDAREA, int IADD, char* out_file);

#endif
