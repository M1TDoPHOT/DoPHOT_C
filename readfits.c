#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "logh.h" 
#include "tuneable.h" 
#include "fitsio.h" 
#include "cast_arr.h" 
#include "readfits.h" 

// read an input 16bit fits file using cfitsio
// cfitsio automatically rescales by bscale and bzero if keywords are set

//inputs: the null terminated name of the fits file to read, 
//pointers for the number of columns nx and rows ny to be updated by cfits io

int** readfits_(char* r_file, int* nx_ptr, int* ny_ptr)
{
     int lverb = tune14_.lverb;
     fitsfile* fptr;
     int fio_err = 0; // for reporting file opening and writing errors
     int file_ver = 0;
     int hdu_num, hdu_type, img_type;
     short int naxis;
     short int nx = 0, ny = 0;
     char NAXIS_comment[FLEN_COMMENT];
     char NAXIS1_comment[FLEN_COMMENT];
     char NAXIS2_comment[FLEN_COMMENT];

     int*  arr1d; 
     int** arr2d;
     long npixels, fpixel[2] = {1, 1};
     int nullval = 0;
     int anynull;
    
     //open the file
     if (lverb > 10){
          fprintf(logfile, 
               "readfits: opening %s for read \n", r_file);
     }
     fits_open_file(&(fptr), r_file, 0, &fio_err);
     if (fio_err != 0){
          printf("readfits_ fitsio error in open_file: %d\n", fio_err);
          fits_close_file(fptr, &fio_err);
          return NULL;
     }

     //get file info to verify file type
     fits_get_num_hdus(fptr, &hdu_num,  &fio_err);
       file_ver += fio_err;
     fits_get_hdu_type(fptr, &hdu_type, &fio_err);
       file_ver += fio_err;
     fits_read_key(fptr, TSHORT, "NAXIS",  &naxis,
                   NAXIS_comment,  &fio_err);
       file_ver += fio_err;
     if (file_ver != 0){
          printf("file verification error: %d\n", file_ver);
          fits_close_file(fptr, &fio_err);
          return NULL;
     }
     if ((hdu_num != 1) || (hdu_type != 0) || (naxis != 2)){
          printf("NOT VALID FILE FOR DOPHOT ANALYSIS\n");
          printf("either too many extensions,\n");
          printf("          invalid HDU type,\n");
          printf("    or NAXIS in header != 2\n");
          return NULL;
     }
     if (lverb > 10){
          fprintf(logfile, "     num HDUs = %d\n", hdu_num );
          fprintf(logfile, "     HDU type = %d\n", hdu_type);
          fprintf(logfile, "     NAXIS    = %d\n", naxis);
     }

     //getting number of axes and data type
     fits_get_img_type(fptr, &img_type, &fio_err);
     fits_read_key(fptr, TSHORT, "NAXIS1",  &nx,
                   NAXIS1_comment,  &fio_err);
     fits_read_key(fptr, TSHORT, "NAXIS2",  &ny,
                   NAXIS2_comment,  &fio_err);
     if ((nx == 0) || (ny == 0)){
          printf("no NAXIS1 or NAXIS2 keywords found\n");
          return NULL;
     }
     if (lverb > 10){
          fprintf(logfile, "     NAXIS1 (nx)  = %d\n", nx);
          fprintf(logfile, "     NAXIS2 (ny)  = %d\n", ny);
          fprintf(logfile, "     data_type    = %d\n", img_type);
     }

     //set the output array and npixels to proper length
     npixels = ny*nx;
     arr1d = malloc_int_1darr(ny* nx);
     arr2d = malloc_int_2darr(ny, nx);

     //read in the image
     fits_read_pix(fptr, TINT, fpixel, npixels, &nullval, 
          arr1d, &anynull, &fio_err);
     if (fio_err != 0){
          printf("readfits_ fitsio error in read_img: %d\n", fio_err);
          fits_close_file(fptr, &fio_err);
          return NULL;
     }

     //recast 1d arr to 2d arr
     recast_int_1dto2darr(ny, nx, arr1d, arr2d);
     free(arr1d);

     fits_close_file(fptr, &fio_err); 
     if (fio_err != 0){
          printf("readfits_ fitsio error in close_file: %d\n", fio_err);
     }          

     *nx_ptr = nx;
     *ny_ptr = ny;

     return arr2d;
}
