#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cast_arr.h"



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
//signed short integer routines

// allocates memory for a 2d c array
signed short int** malloc_ssi_2darr(int size_y, int size_x)
{
     int i, j;
     signed short int** arr ;
     arr = malloc(size_y * sizeof(signed short int*)) ;
     if(arr == NULL){
          fprintf(stderr, "out of memory\n") ;
          return arr;
          }
     for (i = 0; i < size_y; i++){
          arr[i] = malloc(size_x * sizeof(signed short int)) ;
          if(arr[i] == NULL){
               fprintf(stderr, "out of memory\n") ;
               return arr;
               }
          for (j = 0; j < size_x; j++){
               arr[i][j] = 0;
          }
     }
     return arr;
}

// free memory for 2d array
void free_ssi_2darr(int size_y, signed short int **arr)
{
     int i ;
     for(i = 0; i < size_y; i++){
          free(arr[i]);
     }
     free(arr);
}


// converts a 2d fortran array (actually a 1d c array) to a 2d c array 
void recast_ssi_1dto2darr(int ny, int nx, signed short int *arr_ptr, signed short int **arr)
{
     int i,j ;
     for (i = 0; i < ny; i++){
          for (j = 0; j < nx; j++){
               arr[i][j] = *(arr_ptr + nx*i + j) ;
          }
     }
}

// converts a 2d c array to a 2d fortran array which is actually a 1d c array
void recast_ssi_2dto1darr(int ny, int nx, signed short int *arr_ptr, signed short int **arr)
{
     int i,j ;
     for (i = 0; i < ny; i++){
          for (j = 0; j < nx; j++){
               *(arr_ptr + nx*i + j) = arr[i][j];
          }
     }
}

// gets i,j element from a 1d flattened 2d array without an explicit recast
signed short int get_ssi_ij(signed short int* arr_ptr, int nx, int y_elem, int x_elem, int offset)
{
     return *(arr_ptr + nx*y_elem + x_elem + offset) ;
}

// puts i,j element from a 1d flattened 2d array without an explicit recast
void put_ssi_ij(signed short int* arr_ptr, int nx, int y_elem, int x_elem, int offset, signed short int val)
{
     *(arr_ptr + nx*y_elem + x_elem + offset) = val;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
//integer routines

// allocates memory for a 2d c array
int** malloc_int_2darr(int size_y, int size_x)
{
     int i, j;
     int **arr ;
     arr = malloc(size_y * sizeof(int *)) ;
     if(arr == NULL){
          fprintf(stderr, "out of memory\n") ;
          return arr;
          }
     for (i = 0; i < size_y; i++){
          arr[i] = malloc(size_x * sizeof(int)) ;
          if(arr[i] == NULL){
               fprintf(stderr, "out of memory\n") ;
               return arr;
               }
          for (j = 0; j < size_x; j++){
               arr[i][j] = 0;
          }
     }
     return arr;
}

// free memory for 2d array
void free_int_2darr(int size_y, int **arr)
{
     int i ;
     for(i = 0; i < size_y; i++){
          free(arr[i]);
     }
     free(arr);
}


// converts a 2d fortran array (actually a 1d c array) to a 2d c array 
void recast_int_1dto2darr(int ny, int nx, int *arr_ptr, int **arr)
{
     int i,j ;
     for (i = 0; i < ny; i++){
          for (j = 0; j < nx; j++){
               arr[i][j] = *(arr_ptr + nx*i + j) ;
          }
     }
}

// converts a 2d c array to a 2d fortran array which is actually a 1d c array
void recast_int_2dto1darr(int ny, int nx, int *arr_ptr, int **arr)
{
     int i,j ;
     for (i = 0; i < ny; i++){
          for (j = 0; j < nx; j++){
               *(arr_ptr + nx*i + j) = arr[i][j];
          }
     }
}

// gets i,j element from a 1d flattened 2d array without an explicit recast
int get_int_ij(int* arr_ptr, int nx, int y_elem, int x_elem, int offset)
{
     return *(arr_ptr + nx*y_elem + x_elem + offset) ;
}
 
// puts i,j element from a 1d flattened 2d array without an explicit recast
void put_int_ij(int* arr_ptr, int nx, int y_elem, int x_elem, int offset, int val)
{
     *(arr_ptr + nx*y_elem + x_elem + offset) = val;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
// float routines

// allocates memory for a 2d c array
float** malloc_float_2darr(int size_y, int size_x)
{
     int i, j;
     float** arr ;
     arr = malloc(size_y * sizeof(float *)) ;
     if(arr == NULL){
          fprintf(stderr, "out of memory\n") ;
          return arr;
          }
     for (i = 0; i < size_y; i++){
          arr[i] = malloc(size_x * sizeof(float)) ;
          if(arr[i] == NULL){
               fprintf(stderr, "out of memory\n") ;
               return arr;
               }
          for (j = 0; j < size_x; j++){
               arr[i][j] = 0.0f;
          }
     }
     return arr;
}

// free memory for 2d array
void free_float_2darr(int size_y, float **arr)
{
     int i ;
     for(i = 0; i < size_y; i++){
          free(arr[i]);
     }
     free(arr);
}

// converts a 2d fortran array (actually a 1d c array) to a 2d c array 
void recast_float_1dto2darr(int ny, int nx, float *arr_ptr, float **arr)
{
     int i,j ;
     for (i = 0; i < ny; i++){
          for (j = 0; j < nx; j++){
               arr[i][j] = *(arr_ptr + nx*i + j) ;
          }
     }
}

// converts a 2d c array to a 2d fortran array which is actually a 1d c array
void recast_float_2dto1darr(int ny, int nx, float *arr_ptr, float **arr)
{
     int i,j ;
     for (i = 0; i < ny; i++){
          for (j = 0; j < nx; j++){
               *(arr_ptr + nx*i + j) = arr[i][j];
          }
     }
}

// gets i,j element from a 1d flattened 2d array without an explicit recast
float get_float_ij(float* arr_ptr, int nx, int y_elem, int x_elem, int offset)
{
     return *(arr_ptr + nx*y_elem + x_elem + offset) ;
}

// puts i,j element from a 1d flattened 2d array without an explicit recast
void put_float_ij(float* arr_ptr, int nx, int y_elem, int x_elem, int offset, float val)
{
     *(arr_ptr + nx*y_elem + x_elem + offset) = val;
}

// allocates memory for a 3d c array
float*** malloc_float_3darr(int size_z, int size_y, int size_x)
{
     int i;
     float*** arr ;
     arr = malloc(size_z * sizeof(float **)) ;
     if(arr == NULL){
          fprintf(stderr, "out of memory\n") ;
          return arr;
          }
     for (i = 0; i < size_z; i++){
          arr[i] = malloc_float_2darr(size_y, size_x) ;
          if(arr[i] == NULL){
               fprintf(stderr, "out of memory\n") ;
               return arr;
          }
     }
     return arr;
}

     
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
// char routines

// allocates memory for a 2d c array
// each string is 1 mem element longer than chars per str for null char
char** malloc_char_arr(int n_strs, int chars_per_str)
{
     int i ;
     char** arr ;
     arr = malloc(n_strs * sizeof(char *)) ;
     if(arr == NULL){
          fprintf(stderr, "out of memory\n") ;
          return arr;
          }
     for (i = 0; i < n_strs; i++){
          arr[i] = malloc((chars_per_str+1) * sizeof(char)) ;
          if(arr[i] == NULL){
               fprintf(stderr, "out of memory\n") ;
               return arr;
               }
          arr[i][chars_per_str] = '\0' ;
     }
     return arr;
}

// free memory for 2d array
void free_char_arr(int n_strs, char **arr)
{
     int i ;
     for(i = 0; i < n_strs; i++){
          free(arr[i]);
     }
     free(arr);
}

// converts a 2d fortran array (actually a 1d c array) to a 2d c array 
void recast_char_1dto2darr(int n_strs, int chars_per_str, char *arr_ptr, char **arr)
{
     int i,j ;
     for (i = 0; i < n_strs; i++){
          for (j = 0; j < chars_per_str; j++){
               arr[i][j] = *(arr_ptr + chars_per_str*i + j) ;
          }
     }

}

// converts a 2d c array to a 2d fortran array which is actually a 1d c array
void recast_char_2dto1darr(int n_strs, int chars_per_str, char *arr_ptr, char **arr)
{
     int i,j ;
     for (i = 0; i < n_strs; i++){
          for (j = 0; j < chars_per_str; j++){
               *(arr_ptr + chars_per_str*i + j) = arr[i][j];
          }
     }

}

char* initcast_char_2dto1darr(int n_strs, int chars_per_str, char **arr)
{
     char* longstr = malloc(((n_strs*chars_per_str)+1)*sizeof(char));
     memset(longstr, ' ', n_strs*chars_per_str-1);
     int i;
     for (i = 0; i < n_strs; i++){
          if(arr[i][0] != '\0'){
               strncpy((longstr + i*chars_per_str), arr[i], strlen(arr[i]));
          }
          else{
               strncpy((longstr + i*chars_per_str), "  ", 2);
          }
     }
     longstr[n_strs*chars_per_str] = '\0';
     return longstr;
}

// 1darr malloc routines     
signed short int* malloc_ssi_1darr(int size)
{
     int i;
     signed short int* arr ;
     arr = malloc(size * sizeof(signed short int)) ;
     if(arr == NULL){
          fprintf(stderr, "out of memory\n") ;
          return arr;
          }
     for (i = 0; i < size; i++){
          arr[i] = 0;
     }
     return arr;
}

short int* malloc_si_1darr(int size)
{
     int i;
     short int* arr ;
     arr = malloc(size * sizeof(short int)) ;
     if(arr == NULL){
          fprintf(stderr, "out of memory\n") ;
          return arr;
          }
     for (i = 0; i < size; i++){
          arr[i] = 0;
     }
     return arr;
}

int* malloc_int_1darr(int size)
{
     int i;
     int* arr ;
     arr = malloc(size * sizeof(int)) ;
     if(arr == NULL){
          fprintf(stderr, "out of memory\n") ;
          return arr;
          }
     for (i = 0; i < size; i++){
          arr[i] = 0;
     }
     return arr;
}

float* malloc_float_1darr(int size)
{
     int i;
     float* arr ;
     arr = malloc(size * sizeof(float)) ;
     if(arr == NULL){
          fprintf(stderr, "out of memory\n") ;
          return arr;
          }
     for (i = 0; i < size; i++){
          arr[i] = 0.0f;
     }
     return arr;
}

