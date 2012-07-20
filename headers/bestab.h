#ifndef BESTAB_H_
#define BESTAB_H_
 
/* dophot subroutine converted to c void function 3-14-2012 */

/*    The idea is to create a lookuptable for use in "empirical," as
c     opposed to analytic, fits to our subrasters.  In this
c     implementation we use just one star to create the subraster.
c     Ideally one would use all the stars, we'll start using just one,
c     the nth brightest without a saturated pixel within some number of
c     sigma.  We enter this after varipar, so that we can use averages
c     if needed.  One would hope that we've even gone through an improve
c     cycle, but this isn't necessary.  There's a lot of fussing here
c     taking care of a single bad pixel.  Adjacent pairs of bad pixels
c     will kill the template.  Moments of the empircal PSF are computed
c     for use in addlims */

void bestab_( double (*ONESTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr );

#endif
