#ifndef COVARUPD_H_
#define COVARUPD_H_

/* routine to update covariance matrix for object */

void covarupd_(float *fit_covar, float **save_covar, int nfit, int shape_or_improve_call);
     /* shape_or_improve_call = 1 for shape call
                              = 2 for improve call */

#endif
