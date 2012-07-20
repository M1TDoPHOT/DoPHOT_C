#ifndef IMPAPER_H_
#define IMPAPER_H_

/* dophot subroutine converted to c void function 01-29-2012 */
/* 
c
c  This is an independent way of calculating the difference between the
c  aperture and fit mags.  Large values for abs(APPLE[i][2]) for STAR 
c  objects will indicate problems.  This value equals the difference 
c  between the sky and fit mags.
c
*/

void impaper_(int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, int* K_ptr);

#endif
