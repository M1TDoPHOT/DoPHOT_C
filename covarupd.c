#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tuneable.h"
#include "cast_arr.h"
#include "covarupd.h"

/* updates the covariance matrix for a single object from the working 
   'covariance matrix' which contains sqrt variances along the diagonals 
   and normalized covariances on the off diagonal.  updates a saved 
   array of the same. should be called from either shape or improve */

void covarupd_(float *fit_covar, float **save_covar, int nfit, int shape_or_improve_call)
{
     /* shape_or_improve_call = 1 for shape call
                              = 2 for improve call */
    
     int max_zero;
     int i, j;

     /* zero out old save parameters */
     if (shape_or_improve_call == 1){  //shape call
          max_zero = NPMAX;
     }
     else{
          max_zero = tune4_.nfit1;
     }
     for(j = 0; j < max_zero; j++){
          for(i = 0; i < max_zero; i++){
               save_covar[j][i] = 0.0;
          }
     }
     
     /* recasting the 1 dimensional array from the fit into a 2d saved form */
     recast_float_1dto2darr(nfit, nfit, fit_covar, save_covar); 

} 
