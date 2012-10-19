#include <stdlib.h>
#include "logh.h"
#include "tuneable.h"
#include "starlist_struct.h"
#include "search_struct.h"
#include "subraster_struct.h"
#include "crudestat_struct.h"
#include "fitarrays_struct.h"
#include "addobj_struct.h"
#include "chisq.h"
#include "onefit.h"
#include "guess.h"
#include "parinterp.h"
#include "makemask.h"
#include "transmask.h"
#include "offpic.h"
#include "cosmic.h"
#include "oblit.h"
#include "toofaint.h"
#include "toobright.h"
#include "parupd.h"
#include "errupd.h"
#include "addstar.h"
#include "cast_arr.h"
#include "fillerup.h"
#include "isearch.h"

/* dophot int function converted to c int function 03-02-2012 */

int isearch_(double (*ONESTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr)
{

     /* dereference pointers */
     int NFAST = *NFAST_ptr;
     int NSLOW = *NSLOW_ptr;
     /* rename used common blocks */
     short int* IRECT = tune2_.irect;
     short int* KRECT = tune2_.krect;
     int NFIT1        = tune4_.nfit1;
     int NFIT2        = tune4_.nfit2;
     int NIT2         = 2*tune4_.nit;
     int ITOP         = tune3_.itop; //saturated pixel value
     float ENUFF4     = tune10_.enuff4; //if there are enough pixels present
                                        //in the subraster to perform a fit
     int lverb        = tune14_.lverb;
     float* ACC       = tune15_.acc;
     float* ALIM      = tune15_.alim;
     float PIXTHRESH  = tune18_.pixthresh;
     int NTHPIX       = tune18_.nthpix;

     float** STARPAR = starlist_.starpar;
     int* IMTYPE = starlist_.IMTYPE;
     short int** ADDAREA = addobj_.addarea;
     int NSTOT = search_.nstot;
     float THRESH = search_.thresh;
     int NPT = crudestat_.npt;
     /* common blocks only passed to other functions */
     short int** XX = subraster_.xx;
     float* Z  = subraster_.z;
     float* YE = subraster_.ye;
     float*  A = fitarrays_.a;
     float* FA = fitarrays_.fa;
     float* C_ptr = fitarrays_.c; //not recasting as only passed
     float SUM2 = crudestat_.sum2;

     /* substance of funtion begins here */
     short int* JRECT = malloc_si_1darr(2);
     float* ERR   = malloc_float_1darr(NPMAX);
     float* dummy = malloc_float_1darr(NPMAX);
     static int ISUB = -1;

     float PIXT2, BESTSKY, TTHRESH, THRESH2, HIGHSKY, HIGHTHRESH;
     int SKIP, TOOFAINT, COSMIC, TOOBRIGHT, WIPE, HOLE, NUFFPTS;
     int TRANSMASK_SUM2, TRANSMASK_HSKY;
     int I, J;
     float IFLOAT, JFLOAT;
     int NSPREV, ISEARCH;
     int IY, IXMID, ITEST;
     int INOISE, IBIG;
     double parinterp_return, guess1_return;
     float CHI;
     float dum2;
     float DX, DY; //for pass to fxn only
      
     /* FIRST WE GO THROUGH STARSUBTRACTED DATA LOOKING FOR NEW STARS */
     PIXT2   = PIXTHRESH*PIXTHRESH;
     THRESH2 = THRESH*THRESH;
     for(I = 0; I < 2; I++){
          JRECT[I] = IRECT[I];
          if (THRESH > (float)(ITOP/2)){
               IRECT[I] = KRECT[I];
          }
     }
     NSPREV = NSTOT;

     /* update the passmask_.starmask commonblock variable with the star model */
     makemask_(ONESTAR); 

     IXMID = (int)(IXMID/2);
     for(I = 1; I <= NSLOW; I++){ 
          IY = I;
          guess1_return = guess1_(A, dummy, &IXMID, &IY);
          BESTSKY = (float)guess1_return;
          TTHRESH = BESTSKY + THRESH;
          for(J = 1; J <= NFAST; J++){ 
               if ( (J % NTHPIX) == 1 ){
                    if ( J <= (NFAST - (int)(NTHPIX/2)) ){ 
                         JFLOAT = (float)(J + NTHPIX/2); 
                         IFLOAT = (float)(I); 
                         //populate 'dummy' with position specific average star values
                         parinterp_return = parinterp_(&JFLOAT, &IFLOAT, dummy);
                         BESTSKY = (float)parinterp_return;
                    }
                    else{
                         JFLOAT = (float)(J); 
                         IFLOAT = (float)(I); 
                         //populate 'dummy' with position specific average star values
                         parinterp_return = parinterp_(&JFLOAT, &IFLOAT, dummy);
                         BESTSKY = (float)parinterp_return;
                    }
                    TTHRESH = BESTSKY + THRESH;
               } //end J % NTHPIX if
               IBIG   = (int)BIG[I-1][J-1];
               INOISE = NOISE[I-1][J-1];
               if (   ((float)IBIG   >= TTHRESH        ) 
                   && (INOISE        <  MAGIC          )
                   && ((float)INOISE <= (THRESH2/PIXT2)) ){
                    JFLOAT = (float)(J); 
                    IFLOAT = (float)(I); 
                    //populate 'dummy' with position specific average star values
                    parinterp_return = parinterp_(&JFLOAT, &IFLOAT, dummy);
                    HIGHSKY    = (float)parinterp_return;
                    HIGHTHRESH = HIGHSKY + THRESH;
                    TRANSMASK_HSKY = transmask_(BIG, NOISE, &NFAST, &NSLOW,
                                        &J, &I, &HIGHSKY, &dum2);
                    if ( ((float)IBIG >= HIGHTHRESH) && (TRANSMASK_HSKY) ){
                         if (lverb > 20){
                              fprintf(logfile,"Triggering on new Object:%d %d %d\n",
                                               J, I, IBIG);
                         }
                         fillerup_(BIG, NOISE, &J, &I, &NFAST, &NSLOW);
                         //fillerup changes crudestat common vals so reassign
                         NPT  = crudestat_.npt;
                         SUM2 = crudestat_.sum2;
                         ITEST = (int)(ENUFF4*IRECT[0]*IRECT[1]);
                         NUFFPTS = (NPT >= ITEST);
                         if (!NUFFPTS){
                              if (lverb > 20){
                                   fprintf(logfile,"SKIPPING: irect1&2 and NPT = ");
                                   fprintf(logfile,"%d %d %d \n", 
                                           IRECT[0], IRECT[1], NPT);
                              }
                              SKIP = 1;//true
                         }
                         else{
                              TRANSMASK_SUM2 = transmask_(BIG, NOISE, 
                                                      &NFAST, &NSLOW,
                                                      &J, &I, &SUM2, &dum2);
                              if (!TRANSMASK_SUM2){
                                   if (lverb > 20){
                                        fprintf(logfile,"FAILED TRANSMASK ON ");
                                        fprintf(logfile,"AVERAGE SKY\n");
                                   }
                              }
                              else{
                                   guess1_return = guess1_(A, dummy, &J, &I);
                                   CHI = onefit_(ONESTAR, XX, Z, YE, &NPT,
                                                A, FA, C_ptr, &NFIT1, 
                                                ACC, ALIM, &NIT2, 0); 
                                   SKIP = offpic_(A, &J, &I, &NFAST, &NSLOW, 
                                                     &DX, &DY);
                                   if (!SKIP){
                                        NSTOT += 1;
                                        if (lverb > 20){
                                             fprintf(logfile,
                                             "THIS IS STAR NO. %d\n", NSTOT);
                                        }
                                        //update object with 7+ starpar values, 
                                        //even though only 4 fitted.  
                                        //rest obtained from averages in parinterp
                                        parupd_(A, STARPAR[NSTOT-1], J, I, NFIT2);
                                        errupd_(C_ptr, ERR, &NFIT1);
                                        COSMIC    = cosmic_(ONESTAR, 
                                                       BIG, NOISE,
                                                       &NFAST, &NSLOW, 
                                                       STARPAR[NSTOT-1]);
                                        TOOBRIGHT = toobright_( 
                                                       BIG, NOISE,
                                                       &NFAST, &NSLOW, 
                                                       STARPAR[NSTOT-1]);
                                        WIPE = ( (COSMIC) || (TOOBRIGHT) );
                                        if (!WIPE){
                                             TOOFAINT = toofaint_(STARPAR[NSTOT-1], ERR);
                                             if (TOOFAINT){
                                                  IMTYPE[NSTOT-1] = 7;
                                             }
                                        }
                                        if (WIPE){
                                             IMTYPE[NSTOT-1] = 8;
                                             if (COSMIC){ //sub square cosmic rays
                                                  HOLE  = oblit_(ONESTAR, 
                                                          BIG, NOISE,
                                                          &NFAST, &NSLOW, 
                                                          STARPAR[NSTOT-1]);
                                             }
                                             if (TOOBRIGHT){ //sub ellip stars
                                                  HOLE  = oblit_ellipse_(ONESTAR, 
                                                          BIG, NOISE,
                                                          &NFAST, &NSLOW, 
                                                          STARPAR[NSTOT-1]);
                                             }
                                        }
                                        else if (CHI >= 1.0e10f){
                                             SKIP = 1; //true
                                             if (lverb > 20){
                                                  fprintf(logfile,
                                      "FAILED TO CONVERGE: NO ENTRY IN STAR LIST\n");
                                             }
                                             NSTOT -= 1;
                                        }
                                        else{
                                             addstar_(ONESTAR, 
                                                  BIG, NOISE,
                                                  NFAST, NSLOW, 
                                                  STARPAR[NSTOT-1],
                                                  ADDAREA[NSTOT-1],
                                                  ISUB, 
                                                  0, " ", 0, " ");
                                        }
                                   } //end of !SKIP if
                              } //end if transmask_sum2 if/else
                         } //end if !nuffpts if/else
                    }
               }
          } //end J loop
     } //end I loop

     ISEARCH = NSTOT - NSPREV;
     for(I = 0; I < 2; I++){
          IRECT[I] = JRECT[I];
     }

     /* recasting changed pointers and common block vals 
        and freeing alloced mem */
     search_.nstot = NSTOT;

     free(JRECT); 
     free(ERR);
     free(dummy);

     return ISEARCH;

}
                                                  
          
                         
