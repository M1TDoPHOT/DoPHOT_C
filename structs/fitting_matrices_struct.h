#ifndef FITTING_MATRICES_STRUCT_H
#define FITTING_MATRICES_STRUCT_H

// having grown weary of the mallocs
struct{
     float** c_mat; //[npmax, npmax]
     float**    lu; //[npmax, npmax]
     float** b_mat; //[npmax, npmax]
     float* c_list; //[npmax]
     float*      v; //[mmax]
     float*   vsol; //[mmax]
     int*    index; //[mmax]
}fitting_matrices_;

#endif
