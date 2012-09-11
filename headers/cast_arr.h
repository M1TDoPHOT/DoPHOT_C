#ifndef CAST_ARR_H_
#define CAST_ARR_H_

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// for signed short integers

// allocates memory for a 2d c array
signed short int** malloc_ssi_2darr(int size_y, int size_x);

// free memory for 2d array
void free_ssi_2darr(int size_y, signed short int **arr);

// converts a 2d fortran array (actually a 1d c array) to a 2d c array 
void recast_ssi_1dto2darr(int ny, int nx, signed short int *arr_ptr, signed short int **arr);

// converts a 2d c array to a 2d fortran array which is actually a 1d c array
void recast_ssi_2dto1darr(int ny, int nx, signed short int *arr_ptr, signed short int **arr);

// gets i,j element from a 1d flattened 2d array without an explicit recast
signed short int get_ssi_ij(signed short int* arr_ptr, int nx, int y_elem, int x_elem, int offset);

// puts i,j element from a 1d flattened 2d array without an explicit recast
void put_ssi_ij(signed short int* arr_ptr, int nx, int y_elem, int x_elem, int offset, signed short int val);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// for integers

// allocates memory for a 2d c array
int** malloc_int_2darr(int size_y, int size_x);

// free memory for 2d array
void free_int_2darr(int size_y, int **arr);

// converts a 2d fortran array (actually a 1d c array) to a 2d c array 
void recast_int_1dto2darr(int ny, int nx, int *arr_ptr, int **arr);

// converts a 2d c array to a 2d fortran array which is actually a 1d c array
void recast_int_2dto1darr(int ny, int nx, int *arr_ptr, int **arr);

// gets i,j element from a 1d flattened 2d array without an explicit recast
int get_int_ij(int* arr_ptr, int nx, int y_elem, int x_elem, int offset);

// puts i,j element from a 1d flattened 2d array without an explicit recast
void put_int_ij(int* arr_ptr, int nx, int y_elem, int x_elem, int offset, int val);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// for floats

// allocates memory for a 2d c array
float** malloc_float_2darr(int size_y, int size_x);

// free memory for 2d array
void free_float_2darr(int size_y, float **arr);

// converts a 2d fortran array (actually a 1d c array) to a 2d c array 
void recast_float_1dto2darr(int ny, int nx, float *arr_ptr, float **arr);

// converts a 2d c array to a 2d fortran array which is actually a 1d c array
void recast_float_2dto1darr(int ny, int nx, float *arr_ptr, float **arr);

// gets i,j element from a 1d flattened 2d array without an explicit recast
float get_float_ij(float* arr_ptr, int nx, int y_elem, int j_elem, int offset);

// puts i,j element from a 1d flattened 2d array without an explicit recast
void put_float_ij(float* arr_ptr, int nx, int y_elem, int x_elem, int offset, float val);

// allocates memory for a 3d c array
float*** malloc_float_3darr(int size_z, int size_y, int size_x);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// for chars (1d really)

// allocates memory for a 2d c array
char** malloc_char_arr(int n_chars, int size_chars);

// free memory for 2d array
void free_char_arr(int n_chars, char **arr);

// converts a 2d fortran array (actually a 1d c array) to a 2d c array 
void recast_char_1dto2darr(int n_chars, int size_chars, char *arr_ptr, char **arr);

// converts a 2d c array to a 2d fortran array which is actually a 1d c array
void recast_char_2dto1darr(int n_chars, int size_chars, char *arr_ptr, char **arr);
char* initcast_char_2dto1darr(int n_strs, int chars_per_str, char **arr);



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 1d array malloc routines
signed short int* malloc_ssi_1darr(int size);
short int* malloc_si_1darr(int size);
int* malloc_int_1darr(int size);
float* malloc_float_1darr(int size);


#endif
