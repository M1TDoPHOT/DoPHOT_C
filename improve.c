#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "aperlist_struct.h"
#include "estuff_struct.h"
#include "changed_struct.h"
#include "probpass_struct.h"
#include "starlist_struct.h"
#include "search_struct.h"
#include "subraster_struct.h"
#include "crudestat_struct.h"
#include "fitarrays_struct.h"
#include "eoff_struct.h"
#include "unitize_struct.h"
#include "passimp_struct.h"
#include "fixpass_struct.h"
#include "galpass_struct.h"
#include "addobj_struct.h"
#include "model_struct.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "transmask.h"
#include "offpic.h"
#include "toofaint.h"
#include "empiricals.h" //contains oneemp
#include "guess.h"
#include "fillerup.h"
#include "chisq.h"
#include "onefit.h"
#include "elarea.h"
#include "probgal.h"
#include "add_analytic_or_empirical_obj.h"
#include "parupd.h"
#include "errupd.h"
#include "impaper.h"
#include "pgauss.h"
#include "improve.h"

/* dophot int function converted to c int function 03-04-2012 */

int improve_(double (*ONESTAR_7P)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr)
{

     /* dereference pointers */
     int NFAST = *NFAST_ptr;
     int NSLOW = *NSLOW_ptr;
     /* rename used common blocks */
     int*      IMTYPE = starlist_.IMTYPE;
     float**  STARPAR = starlist_.starpar;
     float**    EMPAR = estuff_.empar;
     short int* EMSUB = estuff_.emsub;
     float*     EMERR = estuff_.emerr;
     float**    APPLE = aperlist_.apple;
     float* PROBG = probpass_.probg;
     short int** ADDAREA = addobj_.addarea;
     float XFRAC = eoff_.xfrac;
     float YFRAC = eoff_.yfrac;
     float UFACTOR = unitize_.ufactor;
     int NSTOT   = search_.nstot;
     int I       = passimp_.iii;
//   not renamed because changed frequently and by subroutines
//   crudestat_.npt (int)
//   galpass_.bigfoot (int/ logical)
//   changed_.useold (int/ logical)
//   passimp_.inimp (int/ logical)
//   fixpass_.fixxy (int/ logical)

     short int* IRECT = tune2_.irect;
     short int* KRECT = tune2_.krect;
     float* ARECT = tune2_.arect;
     int NFIT1    = tune4_.nfit1; //usually 4, sky, x, y, intensity
     int NFIT2    = tune4_.nfit2; //usually 7+, full star parameters
     int NIT      = tune4_.nit;
     float EPERDN = tune11_.eperdn;
     float RNOISE = tune11_.rnoise;
     int lverb    = tune14_.lverb;
     float* ACC   = tune15_.acc;
     float* ALIM  = tune15_.alim;
     float APMAGMAXERR = tune19_.apmagmaxerr;
     int FIXPOS   = tune21_.fixpos;
     int EMENAB   = tune22_.emenab;
     int EMPOK    = tune22_.empok;
     int* WHICH_MODEL = model_.which_model;

     /* only passed to probgal and elarea */
     short int** XX = subraster_.xx;
     float* Z  = subraster_.z;
     float* YE = subraster_.ye;
     float*  A = fitarrays_.a;
     float* FA = fitarrays_.fa;
     float* C_ptr = fitarrays_.c;

     //for deciding if addstar should output copies of model images
     // and/ or images Cleaned of neighbors (cf)
     //flags and files [10-13] denote PS, CPS, MPS, and RPS flags
     //and file name roots respectively.  See param_default_c for
     //details
     int w_models = 0, w_cleans = 0;
     int len_mf_names = 0, len_cf_names = 0;
     char **mf_names, **cf_names;
     char *mf_root, *cf_root, *dir;
     int indx_i;

     if ((tune16_.flags[11][0] == 'Y') || (tune16_.flags[11][0] == 'y')){
          dir     = tune16_.files[9];
          cf_root = tune16_.files[11];
          len_cf_names = strlen(cf_root) +strlen(dir) + 10;
          cf_names = malloc_char_arr(NSTOT, len_cf_names);
          for (indx_i = 0; indx_i < NSTOT; indx_i++){
               sprintf(cf_names[indx_i],"%sd%04d%s.fits", dir, indx_i + 1, cf_root);
          }
          w_cleans = 1;
     }
     else{
          cf_names = malloc_char_arr(NSTOT, 1);
          for (indx_i = 0; indx_i < NSTOT; indx_i++){
               cf_names[indx_i][0] = ' ';
          }
     }
     if ((tune16_.flags[12][0] == 'Y') || (tune16_.flags[12][0] == 'y')){
          dir     = tune16_.files[9];
          mf_root = tune16_.files[12];
          len_mf_names = strlen(mf_root) + strlen(dir) + 10;
          mf_names = malloc_char_arr(NSTOT, len_mf_names);
          for (indx_i = 0; indx_i < NSTOT; indx_i++){
               sprintf(mf_names[indx_i],"%sd%04d%s.fits", dir, indx_i + 1, mf_root);
          }
          w_models = 1;
     }
     else{
          mf_names = malloc_char_arr(NSTOT, 1);
          for (indx_i = 0; indx_i < NSTOT; indx_i++){
               mf_names[indx_i][0] = ' ';
          }
     }

     /* substance of function begins here */
     float*  ERR = malloc_float_1darr(NPMAX);    
     float*    B = malloc_float_1darr(NPMAX);    
     short int* JX    = malloc_si_1darr(2);
     short int* JRECT = malloc_si_1darr(2);
     static int IADD = 1;
     static int ISUB = -1;
     static int NFIT0 = 2;
     int SKIP, CONVERGE, SNOK, TRYEM;
     int JMTYPE, NFIT, IIT;
     int IX, IY;
     float DX, DY;
     int logic_ret;
     double chisq_return, onefit_return, elarea_return, oneemp_return;
     float SKY, STARCHI, EMCHI;
     float dum = 0;
     float TOTSTAR, VARSTAR, VARSKY, ERRELEC, ERRDN, SNPRED;
     int K;
     float this_C;
     int zero_dum = 0;
     int which_shape_model;

     double (*ONESTAR)(short int*, float*, float*, int*, int*);

     passimp_.inimp = 1; //true
     TRYEM = (EMENAB && EMPOK);
     for(I = 1; I <= NSTOT; I++){
          K = I-1; //for array indexes

          /* choose which model to add back to the image for
             improve fitting.  Last subtracted by shape or bestab*/
         
          if (WHICH_MODEL[K] == 0){ //normal specified model
               //possibly empirical
               which_shape_model = 0;
               ONESTAR = ONESTAR_7P;
          }
          if (WHICH_MODEL[K] == 1){ //pgaussogaussian
               which_shape_model = 1;
               ONESTAR = &pgauss2d_;
          }

          EMCHI = 2.0e10f;
          fixpass_.fixxy = ( (IMTYPE[K] - (IMTYPE[K] % 10)) == 10);
          fixpass_.fixxy = ( fixpass_.fixxy && FIXPOS ); 
          if (fixpass_.fixxy){
               JMTYPE = IMTYPE[K] - 10;
          }
          else{
               JMTYPE = IMTYPE[K];
          }
          if ( (JMTYPE != 0) && (JMTYPE != 6) && (JMTYPE != 8) ){
               // add back old model
               if (EMSUB[K] == 1) changed_.useold = 1; //true
               add_analytic_or_empirical_obj(ONESTAR,
                     BIG, NOISE, NFAST, NSLOW,
                     STARPAR, ADDAREA, IADD,
                     0, " ", 0, " ", K, (JMTYPE == 2));
               if (EMSUB[K] == 1) changed_.useold = 0; //false

               //temporary override for improve fit and star subtraction
               WHICH_MODEL[K] = 0;
               ONESTAR = ONESTAR_7P;

               /* populate fit matrix A with average star parameters,
                  even for galaxies,
                  don't alter STARPAR */
               SKY = (float)guess3_(A, STARPAR[K], &IX, &IY);
               if (lverb > 20){
                    fprintf(logfile,"IMPROVING STAR # %d AT %d %d \n",
                                     I, IX, IY);
                    fprintf(logfile,"start params are: ");
                    for (indx_i = 0; indx_i < tune4_.nfit2; indx_i++){
                         fprintf(logfile," %10.5f", A[indx_i]);
                    }
                    fprintf(logfile,"\n ");
               }
     
               if (EMSUB[K] == -1){
                    JRECT[0] = IRECT[0];               
                    JRECT[1] = IRECT[1];               
                    IRECT[0] = KRECT[0];               
                    IRECT[1] = KRECT[1]; 
               } 
               // get crude statistics for inten, sky, and x,y val             
               fillerup_(BIG, NOISE, &IX, &IY, &NFAST, &NSLOW); 
               if (EMSUB[K] == -1){
                    IRECT[0] = JRECT[0];               
                    IRECT[1] = JRECT[1]; 
               }              

               if (fixpass_.fixxy){
                    logic_ret = offpic_(A, &IX, &IY,
                                        &NFAST, &NSLOW, 
                                        &DX, &DY);
                    SNOK = !logic_ret;
                    NFIT = NFIT0;
                    IIT  = 2;
               }
               else{
                    /* here I've punted but I could easily 
                       have made a mask from the input
                       subraster */
                    SNOK = transmask_(BIG, NOISE,
                                        &NFAST, &NSLOW, 
                                        &IX, &IY, &SKY, &dum);
                    NFIT = NFIT1;
                    IIT  = NIT;
               }

               if (SNOK){
                    onefit_return = onefit_(ONESTAR, XX, Z, YE,
                                    &crudestat_.npt, A, FA, C_ptr, 
                                    &NFIT, ACC, ALIM, &IIT, WHICH_MODEL[K]);
                    STARCHI = (float)onefit_return;
               }
               else{
                    if (lverb > 20){
                         fprintf(logfile,"snok NG: star #,npt,ix,iy = ");
                         fprintf(logfile,"%d %d %d %d \n",I,crudestat_.npt,IX,IY);
                    }
               }
               logic_ret = offpic_(A, &IX, &IY, &NFAST, &NSLOW, &DX, &DY);
               SKIP = ( (!SNOK) || (logic_ret) );
               if (!SKIP){
                    CONVERGE = (STARCHI < 9.0e9f);
                    if (CONVERGE){
                         if (JMTYPE != 2){
                              //update with most recent average star parameters
                              parupd_(A, STARPAR[K], IX, IY, NFIT2);
                              errupd_(C_ptr, ERR, &NFIT);
                              if (JMTYPE != 3){
                                   logic_ret = toofaint_(STARPAR[K], ERR); 
                                   if (logic_ret){
                                        JMTYPE = 7;
                                   }
                              } //end JMTYPE != 3

                              /* Changed index. */
                              elarea_return = elarea_(A+4, A+5, A+6);
                              TOTSTAR = 2.0f*pi*STARPAR[K][1]
                                        *(float)elarea_return;
                              VARSTAR = TOTSTAR*EPERDN;
                              VARSKY  = ARECT[0]*ARECT[1]
                                        *(SKY*EPERDN + RNOISE*RNOISE);
                              if ( (VARSTAR + VARSKY) > 0.0f ){
                                   ERRELEC = sqrtf(VARSTAR + VARSKY); 
                                   ERRDN   = ERRELEC/EPERDN;
                                   SNPRED  = 1.086f*ERRDN/TOTSTAR;
                                   if (SNPRED < APMAGMAXERR){
                                        /* impaper uses starlist_.starpar 
                                           so need to recast here */
                                        impaper_(BIG, NOISE,
                                                 &NFAST, &NSLOW, &I);
                                   }
                              }
                              /* Changed index. */
                              this_C = get_float_ij(C_ptr, NFIT1, 1, 1, 0);
                              APPLE[K][3] = 1.086f*(1.0f/A[1])*this_C;
                              if (fixpass_.fixxy){
                                   this_C = get_float_ij(C_ptr, NFIT0, 1, 1, 0);
                                   APPLE[K][3] = 1.086f*(1.0f/A[1])*this_C;
                              }
                              /* I'll risk using fa and c again; 
                                 a is needed by probgal */
                              if (TRYEM){
                                   if (lverb > 20){
                                        fprintf(logfile,
                                        "EMPIRICAL FIT TO STAR # %d", I); 
                                        fprintf(logfile,"AT %d %d\n",
                                        IX, IY);
                                   }
                                   if (EMSUB[K] >= 1){
                                        B[0] = EMPAR[K][0]/UFACTOR;
                                        B[1] = EMPAR[K][1]/UFACTOR;
                                        B[2] = EMPAR[K][2] - IX;
                                        B[3] = EMPAR[K][3] - IY;
                                   }
                                   else{
                                        B[0] = A[0];
                                        B[1] = A[1]/10000.0f;
                                        B[2] = A[2] - XFRAC;
                                        B[3] = A[3] - YFRAC;
                                   }
                                   JX[0] = 0;
                                   JX[1] = 0;
                                   oneemp_return = oneemp_(JX, B, FA,
                                                        &NFIT, &zero_dum); 
                                   chisq_return = chisq_(&oneemp_,
                                                XX, Z, YE, &crudestat_.npt,
                                                B, FA, C_ptr, 
                                                &NFIT, ACC, ALIM, &IIT);
                                   EMCHI = (float)chisq_return;
                                   if (EMCHI < 9.0e9f){
                                        EMPAR[K][0] = B[0]*UFACTOR;
                                        EMPAR[K][1] = B[1]*UFACTOR;
                                        EMPAR[K][2] = B[2] + IX;
                                        EMPAR[K][3] = B[3] + IY;
                                        this_C = get_float_ij(C_ptr, NFIT1, 1, 1, 0);
                                        EMERR[K] = 1.086f*(1.0f/B[1])*this_C;
                                        if (fixpass_.fixxy){
                                             this_C = get_float_ij(C_ptr, NFIT0, 1, 1, 0);
                                             EMERR[K] = 1.086f
                                                 *(1.0f/B[1])*this_C;
                                        }
                                   }
                              }// end TRYEM if
                         }// end JMTYPE != 2 if
          
                         /* It doesn't make much sense to ask if it's 
                            bigger than a star if the central intensity 
                            is */
                         if (A[1] > 0.0f){
                              PROBG[K] =  (float)probgal_(ONESTAR,
                                               XX, Z, YE, &crudestat_.npt, A, FA);
                         }
                         else{
                              PROBG[K] = -9999.0f;
                         }
                    }
                    else{
                         if (JMTYPE != 3){
                              JMTYPE = 4;
                         }
                    } // end CONVERGE if/else
                                   
                    if ( (!TRYEM) || (EMCHI > 9.0e9f) 
                         || (EMSUB[K] == -1) ){
                         EMSUB[K] = min(0, EMSUB[K]);
                         if (EMSUB[K] >= 0){
                              EMPAR[K][0] = 0.0f;
                              EMPAR[K][1] = 0.0f;
                              EMPAR[K][2] = 0.0f;
                              EMPAR[K][3] = 0.0f;
                         }
                    }
                    else{
                         EMSUB[K] = 1;
                    }

                    // given new parameters, subtract new model
                    /* for galaxies, need old model type, but 
                       if type 1 or 3, subtract 7+ parameter model
                       regardless of which_model. */
                    if (which_shape_model == 0){ //normal specified model
                         //possibly empirical
                         WHICH_MODEL[K] = 0;
                         ONESTAR = ONESTAR_7P;
                    }
                    if (which_shape_model == 1){ //pgaussogaussian
                         WHICH_MODEL[K] = 1;
                         if ((JMTYPE == 1) || (JMTYPE == 3)){
                              ONESTAR = ONESTAR_7P;
                         }
                         else{
                              ONESTAR = &pgauss2d_;
                         }
                    }

                    add_analytic_or_empirical_obj(ONESTAR, BIG, NOISE,
                                  NFAST, NSLOW, STARPAR, ADDAREA, ISUB,
                                  w_models, mf_names[K], 
                                  w_cleans, cf_names[K],
                                  K, (JMTYPE == 2));

               }
               else{
                    JMTYPE = 6;
                    if (lverb > 20){
                         fprintf(logfile,"DEACTIVATING STAR # %d ", I);
                         fprintf(logfile,"AT %d %d\n", IX, IY);
                         fprintf(logfile,"emchi & emsub = %f %d\n",
                                         EMCHI, EMSUB[K]);
                    }
               } // end !SKIP if/else
          }

          IMTYPE[K] = JMTYPE;
          if (fixpass_.fixxy){
               IMTYPE[K] = JMTYPE + 10;
          }

          //reset which model even if never entered inner loop
          if (which_shape_model == 0){ //normal specified model
               //possibly empirical
               WHICH_MODEL[K] = 0;
          }
          if (which_shape_model == 1){ //pgaussogaussian
               WHICH_MODEL[K] = 1;
          }

          fixpass_.fixxy = 0;//false
          passimp_.iii += 1; //loop this variable as well as cleaner I
     } //end I loop
     passimp_.inimp = 0; //false

     /* recast all changed pointers and common block vars */
     free(ERR);
     free(B);
     free(JX);
     free(JRECT);
     free_char_arr(NSTOT, mf_names);
     free_char_arr(NSTOT, cf_names);

     return 1;
}
