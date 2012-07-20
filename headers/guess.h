#ifndef GUESS1_H_
#define GUESS1_H_

/* dophot function converted to c double function 02-05-2012 */
/*
C:  SERIOUS PROBLEM HERE.  I'D LIKE TO STOP MAKING GUESS1 A FUNCTION.
C:  BUT IT'S CALLED BY ISEARCH, FILLERUP AND DOPHOT.  SO I NEED TO MAKE
C:  CHANGES CONSISTENTLY
*/
/* guess 1 and guess 2 can also be combined, but with a toggle
   they only differ by one line */

double guess1_(float* A, float* STARPAR, int* IX_ptr, int* IY_ptr);
double guess2_(float* A, float* STARPAR, int* IX_ptr, int* IY_ptr);
double guess3_(float* A, float* STARPAR, int* IX_ptr, int* IY_ptr);

#endif
