#include <stdlib.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "free_parking_struct.h"
#include "subraster_struct.h"
#include "starlist_struct.h"
#include "parpred_struct.h"
#include "fitarrays_struct.h"
#include "skyvar_struct.h"
#include "hubvar_struct.h"
#include "model_struct.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "skyfun_plane.h"
#include "skyfun_hub.h"
#include "parinterp.h"
#include "chisq.h"
#include "varipar.h"

/* dophot subroutine converted to c void function 02-24-2012 */
/* if whichmodel = 0, varipar_plane
   if whichmodel = 1, varipar_hub
*/
void varipar_(int* NSTOT_ptr, int* NFAST_ptr, int* NSLOW_ptr, int whichmodel)
{
     /* recasting pointers */
     int NSTOT = *NSTOT_ptr;
     int NFAST = *NFAST_ptr;
     int NSLOW = *NSLOW_ptr;
     /* renaming used common blocks */
     float* SKYPAR;
     if (whichmodel == 0){
          SKYPAR = skyvar_.skypar;
     }
     if (whichmodel == 1){
          SKYPAR = hubvar_.hubpar;
     }
     int lverb  = tune14_.lverb;
     float SKYGUESS = tune1_.skyguess;
     float* AVA = tune15_.ava;
     int MAX_PERF = tune15_.max_perf;
     float* PARMS = parpred_.parms;
     float*  Z = subraster_.z;
     short int** XX = subraster_.xx;
     float* YE = subraster_.ye;
     float*  A = fitarrays_.a;
     float* FA = fitarrays_.fa;
     float* C_ptr = fitarrays_.c; //note, actually 2x2, but not used, only
                                  //passed to chisq, so not recast here
     int*    IMTYPE  = starlist_.IMTYPE;
     float**  STARPAR = starlist_.starpar;
     float**  SHADOW  = starlist_.shadow;
     float**  SHADERR = starlist_.shaderr;
     int*  WHICH_STAR_MODEL = model_.which_model;

     /* substance of subroutine begins here */
     float* PARWT  = free_parking_.npmaxarray_1;
     float* ROOTS  = free_parking_.npmaxarray_2;
     float* WEIGHT = free_parking_.npmaxarray_3;
     float* PARVAL = free_parking_.npmaxarray_4;
     //none of PARWT, ROOTS, WEIGHT, PARVAL are gaurenteed to be 0, 

     static int MINRMS = 5;
     static int ITSKY  = 100;
     static int K[NPMAX-3];
     static int NK = 4;
     static int frst = 1;// true
     // the 7 following are sky model dependent,
     // first 3 for planar
     static int MINSKY = 25;
     static float SKYACC[NPSKY];
     static float SKYLIM[NPSKY];
     // these 4 for hubble
     static int MINHUB = 100;
     static float HACC[NPHUB];
     static float HLIM[NPHUB];
     float* SKYPARINIT;
     // will use only one set depending on model

     int POSITIVE, GOODSTAR, PERFECT, TYPE1, TYPE3, CONV;
     int I, J;
     int NPERF, NGOOD;
     float SKYMAX, XM, YM;
     float TEMP, USQ, SQ, WT;
     double chisq_return, parinterp_return;
     int NSKYFIT_loc; //passed to chisq
     if (whichmodel == 0){
          NSKYFIT_loc = NSKYFIT;
     }
     if (whichmodel == 1){
          NSKYFIT_loc = NHUBFIT; 
          SKYPARINIT = malloc_float_1darr(NPHUB);
     }

     /* initializing arrays */
     if (frst){
          if (tune4_.nfit2 > 7) NK += (tune4_.nfit2 - 7);
          K[0] = 0;
          K[1] = 4;
          K[2] = 5;
          K[3] = 6;
          if (NK > 4){
               for (I = 4; I < NK; I++){
                    K[I] = I + 3;
               }
          }
          
          if (whichmodel == 0){
               for (J = 0; J < NPSKY; J++){
                    SKYACC[J] = 0.001f;
                    SKYLIM[J] = 0.0f;
               }
          }
          if (whichmodel == 1){
               HACC[0] = 0.01f;
               HACC[1] = -0.1f;
               HACC[2] = -0.1f;
               for (J = 3; J < NPHUB; J++){
                    HACC[J] = 0.1f;
               }
               for (J = 0; J < 4; J++){
                    HLIM[J] = -1.0e-5f;
               }
               for (J = 4; J < NPHUB; J++){
                    HLIM[J] = -1.0e-6f;
               }
          }
          frst = 0;//false
     }

     if (whichmodel == 0){
          for(J = 0; J < NPMAX; J++){
               PARMS[J]  = 0.0f;
          }
     }
     if (whichmodel == 1){
          for(J = 0; J < 4; J++){
               PARMS[K[J]]  = 0.0f;
          }
     }

     //none of PARWT, ROOTS, WEIGHT, PARVAL are gaurenteed to be 0, 
     //so initialize all
     for (J = 0; J < NPMAX; J++){
          PARWT[J]  = 0.0f;
          ROOTS[J]  = 0.0f;
          WEIGHT[J] = 0.0f;
          PARVAL[J] = 0.0f;
     }

     /* TWO PASSES, ONE TO COMPUTE EXPECTED VALUE AND ONE TO COMPUTE RMS */
     CONV  = 1; //true
     NPERF = 0;
     NGOOD = 0;
     SKYMAX = -1.0e-10f;
     for(I = 0; I < NSTOT; I++){
          GOODSTAR = (SHADOW[I][0] > 0.0f);
          TYPE1    = ( (IMTYPE[I] == 1) || (IMTYPE[I] == 11) );
          TYPE3    = ( (IMTYPE[I] == 3) || (IMTYPE[I] == 13) );
          PERFECT  = ( (GOODSTAR) && (TYPE1) );
          PERFECT  = ( (PERFECT) && (WHICH_STAR_MODEL[I] == 0));
          GOODSTAR = ( (GOODSTAR) && (TYPE1 || TYPE3) );
          if ( (PERFECT) && (NPERF <= MAX_PERF) ){
               PARWT[0]  += 1.0f;
               PARVAL[0] += STARPAR[I][0];
               for (J = 1; J < NK; J++){
                    TEMP = 1.0f/SHADERR[I][K[J]];
                    PARWT[K[J]] += TEMP;
                    PARVAL[K[J]] += TEMP*SHADOW[I][K[J]];
               }
               NPERF     += 1;
          }
          if (GOODSTAR){
               if (NGOOD < MAXFIL){ //MAXFIL is the dimension of subraster arrays
                    if (STARPAR[I][0] > SKYMAX){
                         SKYMAX = STARPAR[I][0];
                         /* Changed indices. */
                         XM = STARPAR[I][2];
                         YM = STARPAR[I][3];
                    }
                    NGOOD += 1;
                    /* Changed indices. */
                    XX[NGOOD-1][0] = (short int)(STARPAR[I][2] + 0.5f);
                    XX[NGOOD-1][1] = (short int)(STARPAR[I][3] + 0.5f);
                    Z[NGOOD-1]  = STARPAR[I][0];
                    YE[NGOOD-1]  = 1.0f;
               }
          } 
     } //end I loop
     
     if (lverb > 10){
          fprintf(logfile, "# of stars available for computing typical ");
          fprintf(logfile, "SHAPE (Nperf) = %d \n", NPERF);
          fprintf(logfile, "# of stars available for computing ");
          fprintf(logfile, "MODEL SKY (Ngood) = %d \n", NGOOD);
     }

     if (NPERF >= 1){
          for(J = 0; J < NK; J++){
               PARVAL[K[J]] = PARVAL[K[J]]/PARWT[K[J]];
               AVA[K[J]]    = PARVAL[K[J]];
          } 
          if (lverb > 10){
               fprintf(logfile, "WEIGHTED MEANS: ");
               for(J = 0; J < NPMAX; J++){
                    fprintf(logfile, "%f ", PARVAL[J]);
               }
               fprintf(logfile, "\n");
          }

          if (whichmodel == 0){
               if (NGOOD >= MINSKY){
                    if (SKYPAR[0] == 0.0f){
                         SKYPAR[0] = AVA[0];
                         for(J = 1; J < NPSKY; J++){
                              SKYPAR[J] = 0.0f;
                         }
                    }
                    /* Call the appropriate sky function. */
                    /* skyfun_ declared in skyfun_plane.h */
                    chisq_return = chisq_(&skyfun_, XX, Z, YE, &NGOOD, 
                                          SKYPAR, FA, C_ptr, &NSKYFIT_loc, 
                                          SKYACC, SKYLIM, &ITSKY); 
               }
               else{
                    SKYPAR[0] = AVA[0];
                    for(J = 1; J < NPSKY; J++){
                         SKYPAR[J] = 0.0f;
                    }
               } //end NGOOD> MINSKY if/else
          }

          if (whichmodel == 1){
               if (NGOOD >= MINHUB){
                    if (SKYPAR[3] == 0.0f){
                         SKYPAR[0] = SKYGUESS;
                         SKYPAR[1] = XM;
                         SKYPAR[2] = YM;
                         SKYPAR[3] = SKYMAX - SKYPAR[0];
                         SKYPAR[4] = (float)NFAST;
                         SKYPAR[5] = 0.0f;
                         SKYPAR[6] = (float)NSLOW;
                    }
                    for(J = 0; J < NPHUB; J++){
                         SKYPAR[J] = SKYPAR[J];
                    }
                    /* Call the appropriate sky function. */
                    /* hubfun_ declared in skyfun_hub.h */
                    chisq_return = chisq_(&hubfun_, XX, Z, YE, &NGOOD,
                                          SKYPAR, FA, C_ptr, &NSKYFIT_loc,
                                          HACC, HLIM, &ITSKY);
                    CONV = ((float)chisq_return < 1.0e10f);
                    if (!CONV){
                         for(J = 0; J < NPHUB; J++){
                              SKYPAR[J] = SKYPARINIT[J];
                         }
                    }
               }
               else{
                    SKYPAR[0] = SKYGUESS;
                    SKYPAR[1] = (float)NFAST/2.0f;
                    SKYPAR[2] = (float)NSLOW/2.0f;
                    SKYPAR[3] = 0.0f;
                    SKYPAR[4] = (float)NFAST;
                    SKYPAR[5] = 0.0f;
                    SKYPAR[6] = (float)NSLOW;
               }
               if (lverb > 10){
                    fprintf(logfile, "HUBFUN PARAMETERS: \n");
                    for(J = 0; J < NSKYFIT_loc; J++){
                         fprintf(logfile, "%f ", SKYPAR[J]);
                    }
                    fprintf(logfile, "\n");
               }
          }

          if (lverb > 10){
               if (whichmodel == 0){
                    fprintf(logfile, "SKYFUN PARAMETERS: \n");
               }
               if (whichmodel == 1){
                    fprintf(logfile, "HUBFUN PARAMETERS: \n");
               }
               for(J = 0; J < NSKYFIT_loc; J++){
                    fprintf(logfile, "%f ", SKYPAR[J]);
               }
               fprintf(logfile, "\n");
          }
          if (NPERF >= MINRMS){
               for(I = 0; I < NSTOT; I++){
                    PERFECT = ( (IMTYPE[I] == 1) || (IMTYPE[I] == 11) );
                    PERFECT = ( (PERFECT) && (SHADOW[I][0] != 0.0f) );
                    PERFECT = ( (PERFECT) && (WHICH_STAR_MODEL[I] == 0) );
                    if (PERFECT){
                         /* Changed indices. */
                         parinterp_return = parinterp_(STARPAR[I]+2, STARPAR[I]+3, A);
                         for (J = 0; J < 4; J++){
/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:  I HAVE CONVINCED MYSELF THAT THE UNCERTAINTY PER RESIDUAL SQUARED IS
C:  IS 2*DELTA*(O-C).  BUT IF I USE 1/THIS AS A WEIGHT, ACCIDENTAL CANCELLATIONS
C:  GIVE SCREWEY RESULTS.  SO WE'LL TAKE THE LARGER OF O-C**2 AND SHADERR.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
                              USQ = (SHADOW[I][K[J]] - A[K[J]]);
                              SQ  = USQ*USQ;
                              WT  = (1.0f/SHADERR[I][K[J]]) /max1(SHADERR[I][K[J]],SQ);
                              WEIGHT[K[J]] += WT;
                              PARMS[K[J]]  += WT*SQ;
                         }
                    }
               }//end I loop
               POSITIVE = 1; //true
               for(J = 0; J < 4; J++){
                    PARMS[K[J]] = PARMS[K[J]]/WEIGHT[K[J]];
                    POSITIVE = ( POSITIVE && (PARMS[K[J]] >= 0.0f) );
                    if (PARMS[K[J]] >= 0.0f){
                         ROOTS[K[J]] = sqrtf(PARMS[K[J]]);
                    }
                    else{
                         if (lverb > 10){
                              fprintf(logfile, "NEGATIVE SCATTER : \n");
                              fprintf(logfile, "PARAM# & SCATTER = %d %f \n",
                                                K[J], PARMS[K[J]]);
                         } 
                         PARMS[K[J]] = 0.0f;
                         ROOTS[K[J]] = 0.0f;
                    }
               }// end J loop
               if (lverb > 10){
                    fprintf(logfile, "SCATTER:");
                    for(J = 0; J < NPMAX; J++){
                         fprintf(logfile, "%f ", ROOTS[J]);
                    }
                    fprintf(logfile, "\n");
               } 
          }//end NPERF > MINRMS if

     }
     else{
          if (whichmodel == 0){
               SKYPAR[0] = AVA[0];
               for(J = 1; J < NPSKY; J++){
                    SKYPAR[J] = 0.0f;
               }  
          }
          if (whichmodel == 1){
               SKYPAR[0] = SKYGUESS;
               SKYPAR[1] = (float)NFAST/2.0f;
               SKYPAR[2] = (float)NSLOW/2.0f;
               SKYPAR[3] = 0.0f;
               SKYPAR[4] = (float)NFAST;
               SKYPAR[5] = 0.0f;
               SKYPAR[6] = (float)NSLOW;
          }
  
          if (lverb > 10){
               if (whichmodel == 0){
                    fprintf(logfile, "SKYFUN PARAMETERS: \n");
               }
               if (whichmodel == 1){
                    fprintf(logfile, "HUBFUN PARAMETERS: \n");
               }
               for(J = 0; J < NSKYFIT_loc; J++){
                    fprintf(logfile, "%f ", SKYPAR[J]);
               }
               fprintf(logfile, "\n");
          }
     }// end NPERF >=1 if/else

     /* free locally allocated memory */
     if (whichmodel == 1){
          free(SKYPARINIT);
     }
}

void variparplane_(int* NSTOT_ptr, int* NFAST_ptr, int* NSLOW_ptr)
{
     varipar_(NSTOT_ptr, NFAST_ptr, NSLOW_ptr, 0);
}

void variparhub_(int* NSTOT_ptr, int* NFAST_ptr, int* NSLOW_ptr)
{
     varipar_(NSTOT_ptr, NFAST_ptr, NSLOW_ptr, 1);
}

