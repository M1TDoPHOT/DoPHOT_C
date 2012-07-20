#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tuneable.h"
#include "skyvar_struct.h"
#include "hubvar_struct.h"
#include "drfake_struct.h"
#include "median_struct.h"
#include "skyfun_plane.h"
#include "skyfun_hub.h"
#include "cast_arr.h"
#include "parinterp.h"

/* dophot function converted to c double function 02-03-2012 */

double parinterp_(float* X_ptr, float* Y_ptr, float* STARPAR)
{

     /* dereference pointers */
     float X = *X_ptr;
     float Y = *Y_ptr;

     /* rename used common block variables */
     int tNPSKY = NPSKY; //needed so it can be passed by pointer
     int tNPHUB = NPHUB; //needed so it can be passed by pointer
     int NEEDIT = drfake_.needit;
     float* AVA = tune15_.ava   ;
     char** flags = tune16_.flags;
     int** MEDPIX = median_.medpix;

     /* substance of function begins here */
     short int* IXY = malloc_si_1darr(2);
     int* K         = malloc_int_1darr(tune4_.npar - 3);
     float* DUMMY   = malloc_float_1darr(NPMAX);
     int I;
     
     IXY[0] = (short int)(X + 0.5f);
     IXY[1] = (short int)(Y + 0.5f);
     K[0] = 0;
     for(I = 4; I < tune4_.npar; I++){
          K[I-3] = I;
     }

     for (I = 1; I < (tune4_.npar - 3); I++){
          STARPAR[K[I]] = AVA[K[I]];
     }

     NEEDIT = 0; //false
     /* Decide which sky function to use. */
     if (strncmp(flags[1],"PLANE", 5) == 0){
          STARPAR[0] = skyfun_(IXY, skyvar_.skypar, DUMMY, &tNPSKY, &tNPSKY);
     }
     if (strncmp(flags[1],"HUBBL",5) == 0){
          STARPAR[0] = hubfun_(IXY, hubvar_.hubpar, DUMMY, &tNPHUB, &tNPHUB);
     }
     if (strncmp(flags[1],"MEDIA",5) == 0){
          if (IXY[0] < 0){
               IXY[0] = 0;
          }
          if (IXY[1] < 0){
               IXY[1] = 0;
          }
          STARPAR[0] = MEDPIX[IXY[1]][IXY[0]]; 
     }
     NEEDIT = 1; //true

     /* free allocated memory */
     free(IXY)    ;
     free(K)      ;
     free(DUMMY)  ;

     return STARPAR[0];
}


