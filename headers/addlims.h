#ifndef ADDLIMS_H_
#define ADDLIMS_H_

/* dophot subroutine converted to c void function 02-02-2012 */
/* set JRECT size to reflect size of star ellipse (sigma[xy]) */
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   PSEUDOGAUSSIAN EXP(-T**2) = 1/(1 + T**2 + T**4/2 + T**6/6)
  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

void addlims_(float* STARPAR, short int* JRECT);

#endif
