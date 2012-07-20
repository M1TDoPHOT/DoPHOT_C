#ifndef VARIPAR_H_
#define VARIPAR_H_

/* dophot subroutine converted to c void function 02-24-2012 */
/* if whichmodel = 0, varipar_plane
   if whichmodel = 1, varipar_hub
*/
void varipar_(int* NSTOT_ptr, int* NFAST_ptr, int* NSLOW_ptr, int whichmodel);

void variparplane_(int* NSTOT_ptr, int* NFAST_ptr, int* NSLOW_ptr);

void variparhub_(int* NSTOT_ptr, int* NFAST_ptr, int* NSLOW_ptr);

#endif

