#ifndef EMPIRICALS_H_
#define EMPIRICALS_H_

/* dophot function converted to c function 02-07-2012 */

/* twelve point interpolation on both function and x and y derivatives */
/* a(3) and a(4) are the offsets of template (i.e. model) from center
   of current subraster. */
/* xy contains offset of current pixel from center of current raster */
/* note that since a(3) and a(4) change only rarely we can skip a
   lot of stuff. */
/* need to check boundaries of ix0 iy1 et cetera */
/* note that derivatives are now offset by 1/2 */

double oldemp_(short int* xy, float* a, float* fa, int* m_ptr, int* mmax_ptr);

double oneemp_(short int* xy, float* a, float* fa, int* m_ptr, int* mmax_ptr);

#endif
