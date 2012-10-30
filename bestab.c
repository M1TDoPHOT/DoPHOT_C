#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "free_parking_struct.h"
#include "starlist_struct.h"
#include "estuff_struct.h"
#include "search_struct.h"
#include "oimdev_struct.h"
#include "imdev_struct.h"
#include "eoff_struct.h"
#include "empmom_struct.h"
#include "changed_struct.h"
#include "clobber_struct.h"
#include "addobj_struct.h"
#include "model_struct.h"
#include "cast_arr.h"
#include "parinterp.h"
#include "guess.h"
#include "newfits.h"
#include "tagi4.h"
#include "add_analytic_or_empirical_obj.h"
#include "pgauss.h"
#include "bestab.h"
 
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

void bestab_( double (*ONESTAR_7P)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr )
{
     
     /* dereference pointers */
     int NFAST = *NFAST_ptr;
     int NSLOW = *NSLOW_ptr;

     /* rename used common block variables 
        which AREN'T used/changed frequently in subroutines */
     int*      IMTYPE = starlist_.IMTYPE;
     float**  STARPAR = starlist_.starpar; 
//     float**    EMPAR = estuff_.empar;
     short int* EMSUB = estuff_.emsub;
     short int** ADDAREA = addobj_.addarea;
     float CINT  = eoff_.cint;
     float XFRAC = eoff_.xfrac;
     float YFRAC = eoff_.yfrac;
     int NSTOT = search_.nstot;
     int* WHICH_MODEL = model_.which_model;

     float**   DXOMP = oimdev_.dxomp; 
     float**   DYOMP = oimdev_.dyomp;
     int**       OMP = oimdev_.omp; 
     float**   DXEMP = imdev_.dxemp;
     float**   DYEMP = imdev_.dyemp;
     int**       EMP = imdev_.emp; 
//   empmom_.e5 (float)
//   empmom_.e6 (float)
//   empmom_.e7 (float)
//   empmom_.e2 (float)
//   empmom_.e2old (float)
//   changed_.useold (int/logical)
//   clobber_.iold   (int)

     int lverb    = tune14_.lverb;
     char** flags = tune16_.flags;
     char** files = tune16_.files;
     int EMPOK    = tune22_.empok;
     int IHTAB    = tune22_.ihtab;
     float EMINSIG = tune22_.eminsig;
     float TARGX  = tune23_.targx;
     float TARGY  = tune23_.targy;
     float TARGZ  = tune23_.targz;
     int NEMPSKI  = tune23_.nempski;

     /* substance of subroutine begins here */
     int CHOOME, FIRST, IOK;
     int*       IT1 = malloc_int_1darr(NSMAX);
     short int* IT2 = malloc_si_1darr(NSMAX);
     float*       A = free_parking_.npmaxarray_1;

     static int IADD =  1;
     static int ISUB = -1;
//     static int JADD =  2;
//     static int JSUB = -2;
     static int INIT = 1; //true
     static int CLOBBERIOLDNEXT;
     int one_dum = 1; //to pass as ptr to fxn tagi4

     int I, J, K; // K == I-1 for indexing purposes
     int IIND, JIND; // indexes when loop goes from 
                     // -IHSIDE to +IHSIDE for math reasons
     int II, JJ;
     float temp, holder;
     int this_mp1, this_mp2, this_mp3, this_mp4;
     int IBEST, IX, IY;
     int IXL, IYL;
     float SKY, DX, DY, D2;
     int big_iok_flag; //for getting around a fortran goto
     float SUM0, SUMX, SUMY, SUMXX, SUMYY, SUMXY; 
     float AVEX, AVEY, XXMOM, YYMOM, XYMOM, DENOM;
     float T1, T2, T3, T4;
     int ISEE;

     double (*ONESTAR)(short int*, float*, float*, int*, int*);

     for(I = 0; I < NSTOT; I++){
          IT1[I] = 0;
          IT2[I] = 0;
     }

     if (INIT){
          empmom_.e2    = 0.0f;
          empmom_.e2old = 0.0f;
          clobber_.iold = 0;
          CLOBBERIOLDNEXT  = 0;
          INIT          = 0; //false
     }
     else{
          clobber_.iold = CLOBBERIOLDNEXT;
          fprintf(logfile,"old template was %d\n", clobber_.iold);
     }

     IOK = 0;
     EMPOK = 0; //false TUNEABLE CHANGED IN ROUTINE!!!!! FIX ME
     for(I = 1; I <= NSTOT; I++){ //for every object
          K = I-1;
          if (IMTYPE[K] == 1){
               IOK += 1; //IOK counts the stars
               IT1[IOK-1] = (int)(STARPAR[K][1] + 0.5f); //brightness = chance of selection
               IT2[IOK-1] = (short int)(I); //index held by IT2
               if (strncmp(flags[9], "YES", 3) == 0){
                    temp = fabsf(logf(TARGZ/STARPAR[K][1]));
                    CHOOME = (temp <= 1.0f);
                    parinterp_(STARPAR[K]+2, STARPAR[K]+3, A);
                    temp = (TARGX - STARPAR[K][2])*(TARGX - STARPAR[K][2]);
                    CHOOME = ((CHOOME) && (temp < 4.0f*A[4]));
                    temp = (TARGY - STARPAR[K][3])*(TARGY - STARPAR[K][3]);
                    CHOOME = ((CHOOME) && (temp < 4.0f*A[6]));
                    if (CHOOME){
                         IT1[IOK-1] = 2000000000; //selected star has very good chance
                    }
               }
          }
     }// end I loop

     if (IOK >= (1 + NEMPSKI)){ //start with best star by IOK scale
          tagi4_(IT1, &one_dum, &IOK, IT2);
          IOK -= NEMPSKI;
          FIRST = 1; //true
          // move old empirical data to old array holders
          empmom_.e2old = empmom_.e2;
          for (J = -IHSIDE; J <= IHSIDE; J++){
               JIND = J + IHSIDE;
               for (I = -IHSIDE; I <= IHSIDE; I++){
                    IIND = I + IHSIDE;
                    OMP[JIND][IIND] = EMP[JIND][IIND];
                    if ((I != IHSIDE) && (J != IHSIDE)){
                      DXOMP[JIND][IIND] = DXEMP[JIND][IIND];
                      DYOMP[JIND][IIND] = DYEMP[JIND][IIND];
                    }
               }
          }
          if (FIRST){
               FIRST = 0; //false
          }

          big_iok_flag = 1; //start with flag set to true so you enter loop
          while (big_iok_flag){
               big_iok_flag = 0; //unless otherwise toggled, do not loop again
               // assume ok star unless told otherwise
               EMPOK = 1; //true TUNEABLE CHANGED IN ROUTINE!!!!! FIX ME
               IBEST = IT2[IOK-1]; //choose brightest index
               if (WHICH_MODEL[IBEST-1] == 0){ //normal specified model
                    ONESTAR = ONESTAR_7P;
               }
               if (WHICH_MODEL[IBEST-1] == 1){ //pgaussogaussian
                    ONESTAR = &pgauss2d_;
               }

               SKY   = (float)(guess3_(A, STARPAR[IBEST-1], &IX, &IY));
               CINT  = STARPAR[IBEST-1][1];
               XFRAC = STARPAR[IBEST-1][2] - IX;
               YFRAC = STARPAR[IBEST-1][3] - IY;
               // TUNEABLE CHANGED IN ROUTINE!!!!! FIX ME
               // star completely in field?
               EMPOK = ((EMPOK) && ( (IX - IHTAB) >= 1));
               EMPOK = ((EMPOK) && ( (IX + IHTAB) <= NFAST));
               EMPOK = ((EMPOK) && ( (IY - IHTAB) >= 1));
               EMPOK = ((EMPOK) && ( (IY + IHTAB) <= NSLOW));
               // reject if not in field
               if (!EMPOK){
                    IOK -= 1;
                    if (lverb > 1){
                         fprintf(logfile,
                         "star %d too near edge, rejected\n", IBEST);
                    }
                    if (IOK >= 1){
                         big_iok_flag = 1; //true
                         // if IOK is big, reloop and do nothing else in this loop
                         // goes to next best star by IOK standards
                         /* enforced by having !big_iok_flag conditions
                            on all following statements possible entries*/    
                    }
               } 
               else{ // keep star if in field
                    // if formerly an empirical star, add back old empirical model
                    if (EMSUB[IBEST-1] == 1) changed_.useold = 1; //true
                    add_analytic_or_empirical_obj(ONESTAR, 
                          BIG, NOISE, NFAST, NSLOW,
                          STARPAR, ADDAREA, IADD,
                          0, " ", 0, " ", IBEST-1, 0);
                    if (EMSUB[IBEST-1] == 1) changed_.useold = 0; //reset false

                    for (J = -IHTAB; J <= IHTAB; J++){
                         JJ = J + IY - 1;
                         for (I = -IHTAB; I <= IHTAB; I++){
                              II = I + IX - 1;
                              if (NOISE[JJ][II] == MAGIC){
                                   DX = I - A[2];
                                   DY = J - A[3];
                                   D2 = DX*DX/A[4] + 2.0f*A[5]*DX*DY + DY*DY/A[6];
                                   temp = 0.0f;
                                   if (D2 <= EMINSIG*EMINSIG){
                                        EMPOK = 0; //false FIX ME
                                   }
                                   else{
                                        EMPOK = ((EMPOK) && (NOISE[JJ  ][II-1] != MAGIC));
                                        EMPOK = ((EMPOK) && (NOISE[JJ  ][II+1] != MAGIC));
                                        EMPOK = ((EMPOK) && (NOISE[JJ-1][II  ] != MAGIC));
                                        EMPOK = ((EMPOK) && (NOISE[JJ+1][II  ] != MAGIC));

                                        temp += (float)BIG[JJ  ][II-1];
                                        temp += (float)BIG[JJ  ][II+1];
                                        temp += (float)BIG[JJ-1][II  ];
                                        temp += (float)BIG[JJ+1][II  ];
                                        temp = temp/4.0f - SKY;
                                   }
                              }
                              else{
                                        temp = (float)BIG[JJ][II] - SKY;
                              }
                              EMP[J+IHSIDE][I+IHSIDE] = (10000.0f*temp*
                                                   (1.0f/STARPAR[IBEST-1][1]) + 0.5f);
                         }// end I loop
                    }// end J loop

                    if (EMPOK) { // if new empirical template
                                 // force analytic subtraction, even if was
                                 // empirical previously
                         holder = EMSUB[IBEST-1];
                         EMSUB[IBEST-1] = 0; 
                         add_analytic_or_empirical_obj(ONESTAR, 
                               BIG, NOISE, NFAST, NSLOW,
                               STARPAR, ADDAREA, ISUB,
                               0, " ", 0, " ", IBEST-1, 0); 
                         EMSUB[IBEST-1] = holder;
                    }
                    else{ 
                         if (EMSUB[IBEST-1] >= 1) changed_.useold = 1; //true
                         add_analytic_or_empirical_obj(ONESTAR, 
                               BIG, NOISE, NFAST, NSLOW,
                               STARPAR, ADDAREA, ISUB,
                               0, " ", 0, " ", IBEST-1, 0); 
                         if (EMSUB[IBEST-1] >= 1) changed_.useold = 0; //false
                    } 

                    if (!EMPOK){ 
                         IOK -= 1;
                         if (lverb > 1){
                              fprintf(logfile,
                              "star %d rejected as template\n", IBEST);
                         }
                         if (IOK >= 1){
                              big_iok_flag = 1; //true, will reloop 
                              // using next best star
                         }
                    }
                    else{
                         IXL = (int)(2.0f*sqrtf(STARPAR[IBEST-1][4])) + 1;
                         IYL = (int)(2.0f*sqrtf(STARPAR[IBEST-1][6])) + 1;
                         SUM0 = 0.0f;
                         SUMX = 0.0f;
                         SUMY = 0.0f;
                         SUMXX = 0.0f;
                         SUMYY = 0.0f;
                         SUMXY = 0.0f;
                         for (J = -IYL; J <= IYL; J++){
                              JIND = J + IHSIDE;
                              for (I = -IXL; I <= IXL; I++){
                                   IIND = I + IHSIDE;
                                   SUM0  += (float)((EMP[JIND][IIND]));
                                   SUMX  += (float)((EMP[JIND][IIND])*I); 
                                   SUMY  += (float)((EMP[JIND][IIND])*J); 
                                   SUMXX += (float)((EMP[JIND][IIND])*I*I); 
                                   SUMXY += (float)((EMP[JIND][IIND])*I*J); 
                                   SUMYY += (float)((EMP[JIND][IIND])*J*J); 
                              }
                         }
                         AVEX = SUMX/SUM0;
                         AVEY = SUMY/SUM0;
                         XXMOM = SUMXX/SUM0 - AVEX*AVEX;
                         YYMOM = SUMYY/SUM0 - AVEY*AVEY;
                         XYMOM = SUMXY/SUM0 - AVEX*AVEY;
                         DENOM = XXMOM*YYMOM - XYMOM*XYMOM;
                         empmom_.e5 =  YYMOM/DENOM;
                         empmom_.e7 =  XXMOM/DENOM;
                         empmom_.e6 = -XYMOM/DENOM;
                         empmom_.e5 = 1.0f/empmom_.e5; 
                         empmom_.e7 = 1.0f/empmom_.e7; 
                         for (J = -IHSIDE; J < IHSIDE; J++){
                              JIND = J + IHSIDE;
                              for (I = -IHSIDE; I < IHSIDE; I++){
                                   IIND = I + IHSIDE;
                                   this_mp1 = EMP[JIND+1][IIND+1];
                                   this_mp2 = EMP[JIND+1][IIND  ];
                                   this_mp3 = EMP[JIND  ][IIND+1];
                                   this_mp4 = EMP[JIND  ][IIND  ];
                                   T1 = (float)(this_mp1 - this_mp2);
                                   T2 = (float)(this_mp3 - this_mp4);
                                   T3 = (float)(this_mp1 - this_mp3);
                                   T4 = (float)(this_mp2 - this_mp4);
                                   DXEMP[JIND][IIND] = (T1 + T2)/2.0f;
                                   DYEMP[JIND][IIND] = (T3 + T4)/2.0f;
                              }
                         }

                         if (lverb > 1){
                              fprintf(logfile,
                                      "star %d used as template\n", IBEST);
                              fprintf(logfile,"avex & avey = %f %f \n",
                                      AVEX, AVEY); 
                              fprintf(logfile,"e(5-7) = %f %f %f\n",
                                      empmom_.e5, empmom_.e6, empmom_.e7);
                         }

                         EMSUB[IBEST-1] = -1; // empirical template star now flagged
                         empmom_.e2 = STARPAR[IBEST-1][1];

                         /* if old template was a different star, add the analytic psf 
                            previously subtracted back to the image, and subtract 
                            the old empirical template (itself) in prep for improve.
                            deactivate as empirical psf star */
                         /* additionally, if EMPIRICAL template is new, flag with 0
                            so old template is added back on next pass, not analytic */
                         if ((clobber_.iold != 0) && (clobber_.iold != IBEST)){
                              fprintf(logfile,"deactivating old empirical psf \n");
                              fprintf(logfile,"old template was %d\n", 
                                               clobber_.iold);
                   
                              // add back old analytic
                              add_analytic_or_empirical_obj(ONESTAR, 
                                   BIG, NOISE, NFAST, NSLOW,
                                   STARPAR, ADDAREA, 1,
                                   0, " ", 0, " ", clobber_.iold-1, 0);
                              EMSUB[clobber_.iold-1] = 1; // empirical template unflagged

                              // subtract out old empirical, should be perfect
                              changed_.useold = 1; //true
                              add_analytic_or_empirical_obj(ONESTAR, 
                                   BIG, NOISE, NFAST, NSLOW,
                                   STARPAR, ADDAREA, -1,
                                   0, " ", 0, " ", clobber_.iold-1, 0);
                              changed_.useold = 0; //false

                         } 

                         CLOBBERIOLDNEXT = IBEST;
// HUGE BUG now corrected                       clobber_.iold = IBEST;

                         if (strncmp(flags[8], "YES", 3) == 0){
                              ISEE = (2*IHSIDE) + 1;
                              newfits_(ISEE, ISEE, EMP, files[7], 0, " ");
                         }
                    } //end !EMPOK if/else
               } // other end !EMPOK if/else
          } //end while big IOK which was the goto get around
     } //end IOK > 1+NEMPSKI if

     if (!EMPOK){
          if (lverb > 1){
               fprintf(logfile,"*** no suitable templates ***\n");
          }
     }

     /* restore changed common blocks and free locally allocated mem */
     eoff_.cint  = CINT;
     eoff_.xfrac = XFRAC;
     eoff_.yfrac = YFRAC;
     tune22_.empok = EMPOK; //FIX ME tuneable changed

     free(IT1);
     free(IT2);

}
                               
               

