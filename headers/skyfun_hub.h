#ifndef _SKYFUN_HUB_H_
#define _SKYFUN_HUB_H_

/* dophot function converted to c double function 02-03-2012 */
/*
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:  NOTE THAT NONE OF THE QUANTITIES ARE SCALED.  IN PRINCIPAL WE COULD
C:  SCALE BOTH THE INTENSITIES AND THE POSITIONS, THOUGH THE LATTER LOOKS
C:  TRICKY TO ME.  NOTE THAT DERIVATIVES OF 5, 6 AND 7 ARE COMPUTED BUT
C:  NOT FIT FOR AT THE MOMENT.  WE ADOPT THE SHAPE OF THE SEEING PROFILE
C:  BUT COULD EASILY SOLVE FOR ALL 7 PARAMETERS IN VARIPAR.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/

double hubfun_(short int* IX, float* A, float* FA, int* M_ptr, int* MMAX_ptr);

#endif
