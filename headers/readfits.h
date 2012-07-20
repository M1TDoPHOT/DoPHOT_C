#ifndef READFITS_H_
#define READFITS_H_

// read an input 16bit fits file using cfitsio
// cfitsio automatically rescales by bscale and bzero if keywords are set

//inputs: the null terminated name of the fits file to read, 
//pointers for the number of columns nx and rows ny to be updated by cfits io

int** readfits_(char* r_file, int* nx_ptr, int* ny_ptr);

#endif
