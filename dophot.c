#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "logh.h"
#include "tuneable.h"
#include "starlist_struct.h"
#include "aperlist_struct.h"
#include "unitize_struct.h"
#include "undergnd_struct.h"
#include "ctimes_struct.h"
#include "drfake_struct.h"
#include "search_struct.h"
#include "probpass_struct.h"
#include "median_struct.h"
#include "trans7_struct.h"
#include "fixpass_struct.h"
#include "estuff_struct.h"
#include "subraster_struct.h"
#include "byvirtue_struct.h"
#include "fitarrays_struct.h"
#include "fitting_matrices_struct.h"
#include "free_parking_struct.h"
#include "hubvar_struct.h"
#include "imdev_struct.h"
#include "oimdev_struct.h"
#include "parpred_struct.h"
#include "passmask_struct.h"
#include "skyvar_struct.h"
#include "addobj_struct.h"
#include "model_struct.h"
#include "cast_arr.h"
#include "tuneup.h"
#include "readfits.h"
#include "makenoise.h"
#include "warmstart.h"
#include "pgauss.h"
#include "gauss.h"
//#include "galtmodel.h"
//#include "pgaltmodel.h"
//#include "sersic.h"
#include "extpgauss.h"
#include "medfil.h"
#include "varipar.h"
#include "makemask.h"
#include "shape.h"
#include "improve.h"
#include "isearch.h"
#include "paravg.h"
#include "io_dophotc.h"
#include "outputs.h" //contains stdotpt and badotpt
#include "bestab.h"
#include "write_to_log.h"
#include "newfits.h"
#include "stampout.h"

// input is optionally the parameter modifications file
// if it is not given, it will be prompted for 
// in paramfile.c via tuneup.c
int main( int argc, char* argv[])
{
     clock_t starttime, endtime;
     starttime = clock();
     printf("Welcome to DoPHOT v 4.2\n");

     /* allocating memory for the non-tuneable common blocks */
     float** STARPAR   = malloc_float_2darr(NSMAX, NPMAX);
     starlist_.starpar = STARPAR;
     float** SHADOW    = malloc_float_2darr(NSMAX, NPMAX);
     starlist_.shadow  = SHADOW;
     float** SHADERR   = malloc_float_2darr(NSMAX, NPMAX);
     starlist_.shaderr = SHADERR;
     float*** SHADCOVAR  = malloc_float_3darr(NSMAX, NPMAX, NPMAX);
     starlist_.shadcovar = SHADCOVAR;
     int* IMTYPE       = malloc_int_1darr(NSMAX);
     starlist_.IMTYPE  = IMTYPE;

     float** APPLE     = malloc_float_2darr(NSMAX, NAPPLE);
     aperlist_.apple   = APPLE;
     
     short int** ADDAREA = malloc_ssi_2darr(NSMAX, 4);
     addobj_.addarea   = ADDAREA;

     float* PROBG      = malloc_float_1darr(NSMAX);
     probpass_.probg   = PROBG;

     int** MEDPIX; //only malloced if MEDIAN sky

     short int* EMSUB  = malloc_ssi_1darr(NSMAX);
     estuff_.emsub     = EMSUB;
     float** EMPAR     = malloc_float_2darr(NSMAX, NPMAX); 
                         //NPMAX, not NEP because must be compatible with
                         //parinterp and guess functions
     estuff_.empar     = EMPAR;
     float* EMERR      = malloc_float_1darr(NSMAX);
     estuff_.emerr     = EMERR;

     float* ZE         = malloc_float_1darr(MAXFIL);
     subraster_.z      = ZE;
     short int** XX    = malloc_ssi_2darr(MAXFIL, 2);
     subraster_.xx     = XX;
     float* YE         = malloc_float_1darr(MAXFIL);
     subraster_.ye     = YE;

     float* chi_       = malloc_float_1darr(4);
     byvirtue_.chi     = chi_;

     int* which_model  = malloc_int_1darr(NSMAX);
     model_.which_model= which_model;
     int* tested       = malloc_int_1darr(NSMAX);
     model_.tested     = tested;

     float*  a         = malloc_float_1darr(NPMAX);
     fitarrays_.a      = a;
     float* fa         = malloc_float_1darr(NPMAX);
     fitarrays_.fa     = fa;
     float*  c         = malloc_float_1darr((2*NPMAX)*(2*NPMAX));
     fitarrays_.c      = c;
     float*  b         = malloc_float_1darr(2*NPMAX);
     fitarrays_.b      = b;
     float* fb         = malloc_float_1darr(2*NPMAX);
     fitarrays_.fb     = fb;

     float** c_mat           = malloc_float_2darr(NPMAX+1, NPMAX+1);
     fitting_matrices_.c_mat = c_mat;
     float** b_mat           = malloc_float_2darr(NPMAX+1, NPMAX+1);
     fitting_matrices_.b_mat = b_mat;
     float** lu              = malloc_float_2darr(NPMAX+1, NPMAX+1);
     fitting_matrices_.lu    = lu;
     float*  c_list          = malloc_float_1darr(NPMAX);
     fitting_matrices_.c_list = c_list;
     float*  v               = malloc_float_1darr(MMAX);
     fitting_matrices_.v     = v;
     float* vsol             = malloc_float_1darr(MMAX);
     fitting_matrices_.vsol  = vsol;
     int* index_list         = malloc_int_1darr(MMAX);
     fitting_matrices_.index = index_list;

     float*        npmaxarray_1 = malloc_float_1darr(NPMAX);
     free_parking_.npmaxarray_1 = npmaxarray_1;
     float*        npmaxarray_2 = malloc_float_1darr(NPMAX);
     free_parking_.npmaxarray_2 = npmaxarray_2;
     float*        npmaxarray_3 = malloc_float_1darr(NPMAX);
     free_parking_.npmaxarray_3 = npmaxarray_3;
     float*        npmaxarray_4 = malloc_float_1darr(NPMAX);
     free_parking_.npmaxarray_4 = npmaxarray_4;
     short int*    sifourarray_1 = malloc_si_1darr(4);
     free_parking_.sifourarray_1 = sifourarray_1;
     short int*    sifourarray_2 = malloc_si_1darr(4);
     free_parking_.sifourarray_2 = sifourarray_2;
     int*          intfourarray_1 = malloc_int_1darr(4);
     free_parking_.intfourarray_1 = intfourarray_1;
     int*          intfourarray_2 = malloc_int_1darr(4);
     free_parking_.intfourarray_2 = intfourarray_2;
     float**       npmaxbynpmaxarray = malloc_float_2darr(NPMAX, NPMAX);
     free_parking_.npmaxbynpmaxarray = npmaxbynpmaxarray;

     float* hubpar     = malloc_float_1darr(NPHUB);
     hubvar_.hubpar    = hubpar;

     float** dxemp     = malloc_float_2darr((2*IHSIDE),(2*IHSIDE));
     imdev_.dxemp      = dxemp;
     float** dyemp     = malloc_float_2darr((2*IHSIDE),(2*IHSIDE));
     imdev_.dyemp      = dyemp;
     int** emp         = malloc_int_2darr(((2*IHSIDE)+1),((2*IHSIDE)+1));
     imdev_.emp        = emp;

     float** dxomp     = malloc_float_2darr((2*IHSIDE),(2*IHSIDE));
     oimdev_.dxomp     = dxomp;
     float** dyomp     = malloc_float_2darr((2*IHSIDE),(2*IHSIDE));
     oimdev_.dyomp     = dyomp;
     int** omp         = malloc_int_2darr(((2*IHSIDE)+1),((2*IHSIDE)+1));
     oimdev_.omp       = omp;

     float* parms      = malloc_float_1darr(NPMAX);
     parpred_.parms    = parms;

     float** starmask  = malloc_float_2darr(17,17);
     passmask_.starmask = starmask;

     float* skypar     = malloc_float_1darr(NPSKY);
     skyvar_.skypar    = skypar;
     
     /* listing other commonblock vars used used here, 
        and their types, for convenience */
     //unitize_.ufactor (float)
     //undergnd_.NFAST NSLOW; (int)
     //ctimes_.chiimp apertime, filltime, addtime (float)
     //drfake_.needit (int)
     //search_.nstot (int)
     //search_.thresh (float)
     //trans7_.test7 (int)
     //fixpass_.fixxy (int)

     /* substance begins here */
     /* initializations */
     int first, WARM;
     int PASSCNT, NPPP;

     drfake_.needit   = 1; //true     
     fixpass_.fixxy   = 0; //false 
     trans7_.test7    = 0;
     unitize_.ufactor = 100.0f;
     search_.nstot    = 0;    

     /* initialize all of the tuneable.h variables
        by reading in from modified parameter file and from 
        there the default parameter file */

     int command_line_pm;
     char* pm_file_name;
     if (argc != 2){ /* parameter modification file to be 
                        entered manually */
          command_line_pm = 0;
          pm_file_name = " ";
     }
     else{
          command_line_pm = 1;
          pm_file_name = argv[1];
     }
     tuneup_(command_line_pm, pm_file_name); 

     /*   useful recasts: */
     int lverb = tune14_.lverb;
     char** files = tune16_.files;
     char** flags = tune16_.flags;

     //chose model based on flag
     double (*model2d)(short int*, float*, float*, int*, int*);
     double (*model4d)(short int*, float*, float*, int*, int*);
     if (strncmp(flags[0], "PGAUSS", 5) == 0){
          model2d = &pgauss2d_;
          model4d = &pgauss4d_;
     }
     if (strncmp(flags[0], "GAUSS",  5) == 0){
          model2d = &gauss2d_;
          model4d = &gauss4d_;
     }
//     if (strncmp(flags[0], "PGALTMODEL",  5) == 0){
//          model2d = &pgaltmodel2d_;
//          model4d = &pgaltmodel4d_;
//     }
//     if (strncmp(flags[0], "GALTMODEL",  5) == 0){
//          model2d = &galtmodel2d_;
//          model4d = &galtmodel4d_;
//     }
     if (strncmp(flags[0], "SERSIC", 5) == 0){
          printf("SERSIC MODEL CURRENTLY NOT FUNCTIONAL\n");
          printf("USING GAUSSIAN MODEL INSTEAD\n");
          // need seperate upper and lower limits to hold
          // sersic index above zero
          model2d = &gauss2d_;
          model4d = &gauss4d_;
          tune4_.npar  -= 1;
          tune4_.nfit2 -= 1;
          free(flags[0]);
          flags[0] = malloc((strlen("GAUSS") + 1) * sizeof(char));
          strcpy(flags[0], "GAUSS");
//          model2d = &sersic2d_;
//          model4d = &sersic4d_;
     }
     if (strncmp(flags[0], "EXTPGAUSS",  5) == 0){
          model2d = &extpgauss2d_;
          model4d = &extpgauss4d_;
     }

     /* Open log file if desired. */
     char* logname;
     if ((lverb > 0) && (strncmp(files[5],"TERM",4) != 0)){ 
          logname = stringstrip_(files[5]);
          openlog_(logname);
          free(logname);
     }
     // the logfile is now open and accessable 
     // to all files which #include logh.h
     lverb = lverb*10 + 1;
     tune14_.lverb = lverb;
           
     /* Open input picture. */
     int nx, ny;
     int** big = readfits_(files[0], &nx, &ny);          
     int nfast = undergnd_.NFAST = nx; 
     int nslow = undergnd_.NSLOW = ny; 
     
     if (tune18_.nthpix <= 0); tune18_.nthpix = (int)sqrtf((float)nx);
    
     if (strncmp(flags[1], "MEDIA", 5) == 0){
          if(lverb >= 10){ 
               fprintf(logfile, "Median box half-size (X,Y): %d %d\n",
                                tune20_.jhxwid, tune20_.jhywid);
               fprintf(logfile, "Median precision in DN : %d\n",
                                tune20_.mprec);
          }
          PASSCNT = 0;
          MEDPIX = malloc_int_2darr(nslow, nfast);
          median_.medpix    = MEDPIX;
     }

     if(lverb >= 10){ 
          fprintf(logfile, "NTHPIX (initial sky update freq. in pixels) : %d\n",
                           tune18_.nthpix);
     }

     /* making noise array and starting the warmstart if flagged */
     int** noise = malloc_int_2darr(nslow, nfast);
     makenoise_(big, noise, &nfast, &nslow);
     WARM = (flags[5][0] == 'Y');
     char* warmstartfile;
     if (WARM){
          warmstartfile = stringstrip_(files[2]);
          //printf("here before warmstart\n");
          warmstart_(model2d, big, noise, &nfast, &nslow, warmstartfile);
          //printf("here after warmstart\n");
          free(warmstartfile);
     }

     search_.thresh = tune3_.tmin;
     NPPP = 1;
     while (search_.thresh*powf(2, tune3_.tfac) <= tune3_.tmax){
          NPPP += 1;
          search_.thresh = search_.thresh*powf(2, tune3_.tfac);
     }

     int outfiletype; //INTERNAL = 0;  COMPLETE = 1;  INCOMPLETE = 2
     if (strncmp(flags[2], "INTER", 5) == 0) outfiletype = 0;
     if (strncmp(flags[2], "COMPL", 5) == 0) outfiletype = 1;
     if (strncmp(flags[2], "INCOM", 5) == 0) outfiletype = 2;

     int nrmax_dum = 0; //0
     int ncmax_dum = 0; //0
     int napple_dum = NAPPLE;
     int NSTAR;
     int I, K, J;
     char* outputline;
     char* rwa = "w";
     int ierr;
     char* outfilename;
     FILE* outfile;
     char* errfilename; //only used if outfiletype = 0
     FILE* errfile; //only used if outfiletype = 0
     char* shadowfilename;
     FILE* shadowfile;
     char* covarfilename;
     FILE* covarfile;
    
      
     first = ( (tune21_.fixpos) || (WARM) );
     int last = 0;
     while ( ((search_.thresh / tune3_.tmin) >= 0.999f) || (last) ){

          tune22_.emenab = (search_.thresh <= tune22_.emthrsh);
          if (!last){
               if (lverb > 10){
                    fprintf(logfile," \n");        
                    fprintf(logfile,"Starting loop at threshold level %f \n",
                                     search_.thresh);
               }
          }
          else{
               if (lverb > 10){
                    fprintf(logfile," \n");        
                    fprintf(logfile,"Last Loop \n");
               }
          }
          
 
//        printf("here before sky\n");
          if (strncmp(flags[1], "MEDIA", 5) == 0){
               if ((PASSCNT == 0) || (PASSCNT == (NPPP-2))){
                    if (lverb >= 10){
                         fprintf(logfile,"(re-)Making median background picture \n");
                    }
                    medfil_(&nfast, &nslow, big, MEDPIX, 
                            &tune3_.ibot, &tune3_.itop, 
                            &tune20_.jhxwid, &tune20_.jhywid,
                            &nrmax_dum, &ncmax_dum, &tune20_.mprec);
               }
               else{
                    if (lverb >= 10){
                         fprintf(logfile,"Skipping median filtering on this pass \n");
                    }
               }
               PASSCNT += 1;
               variparplane_(&search_.nstot, &nfast, &nslow);
          }
          if (strncmp(flags[1], "PLANE", 5) == 0){
               variparplane_(&search_.nstot, &nfast, &nslow);
          }
          if (strncmp(flags[1], "HUBBL", 5) == 0){
               variparhub_(&search_.nstot, &nfast, &nslow);
          }
//        printf("here after sky\n");
               
          /* need varipar because improve calls guess3 which calls skyfun
             need makemask because improve calls snok which needs mask
             empirical PSF won't be fitted on this pass because emenab
             is false */

//        printf("here before makemask\n");
          makemask_(model2d);
//        printf("here before makemask\n");

          if ((!first) && (!last)){
//             printf("here before isearch\n");
               NSTAR = isearch_(model2d, big, noise, &nfast, &nslow);
//             printf("here after isearch\n");
          }

//        printf("here before shape\n");
          shape_(model2d, model4d, big, noise, &nfast, &nslow);
//        printf("here after shape\n");
//        printf("here before paravg\n");
          paravg_();
//        printf("here after paravg\n");

          if (search_.nstot >= 1){
               if (strncmp(flags[1], "MEDIA", 5) == 0){
                    variparplane_(&search_.nstot, &nfast, &nslow);
               }
               if (strncmp(flags[1], "PLANE", 5) == 0){
                    variparplane_(&search_.nstot, &nfast, &nslow);
               }
               if (strncmp(flags[1], "HUBBL", 5) == 0){
                    variparhub_(&search_.nstot, &nfast, &nslow);
               }
          }

          // improving for empiricals         
          if (tune22_.emenab){
//        printf("here before bestab\n");
                bestab_(model2d, big, noise, &nfast, &nslow);
//        printf("here after bestab\n");
          }
//        printf("here before improve\n");
          improve_(model2d, big, noise, &nfast, &nslow);
//        printf("here after improve\n");

          if ((!first) && (!last)){
               if (lverb > 10){
                    fprintf(logfile,"Ending loop at threshold level %f \n",
                                     search_.thresh);
                    fprintf(logfile,
                    "Number of new objects found on this threshold = %d \n",
                                NSTAR);        
                    fprintf(logfile,
                    "Total number of objects found so far = %d \n", 
                                     search_.nstot);
               }
          }

          if ((!first) && (!last)){
               search_.thresh = search_.thresh/powf(2, tune3_.tfac);
          }

          if ( ((search_.thresh / tune3_.tmin) < 0.999f) && (!last)){
               if (tune3_.dofinalfit){ 
                    last = 1;
               }
          }
          else{
               last = 0; //only last once if ever
          }

          /* outputing results to files */
          outfilename = stringstrip_(files[3]);
          outfile = openac_(&ierr, outfilename, rwa);
          switch(outfiletype){
               case 0:
                    for (I = 0; I < search_.nstot; I++){
                         K = I+1;
                         outputline = sumout_(&K, IMTYPE+I, STARPAR[I], 
                                &tune4_.npar, APPLE[I], &napple_dum, 
                                PROBG+I, which_model[I]);
                         fputs(outputline, outfile);
                    }
                    break;
               case 1:
                    for (I = 0; I < search_.nstot; I++){
                         K = I+1;
                         outputline = stdotpt_(&K, IMTYPE+I, STARPAR[I], 
                                &tune4_.npar, APPLE[I], &napple_dum, 
                                PROBG+I, which_model[I]);
                         fputs(outputline, outfile);
                    }
                    break;
               case 2:
                    for (I = 0; I < search_.nstot; I++){
                         K = I+1;
                         outputline = badotpt_(&K, IMTYPE+I, STARPAR[I], 
                                &tune4_.npar, APPLE[I], &napple_dum,
                                PROBG+I);
                         fputs(outputline, outfile);
                    }
                    break;
               default:
                    for (I = 0; I < search_.nstot; I++){
                         K = I+1;
                         outputline = sumout_(&K, IMTYPE+I, STARPAR[I], 
                                &tune4_.npar, APPLE[I], &napple_dum,
                                PROBG+I, which_model[I]);
                         fputs(outputline, outfile);
                    }
                    break;
          }//end switch
          fclose(outfile);
          free(outfilename);

          /* write shadow image if flagged */
          if (flags[3][0] == 'Y'){
               shadowfilename = stringstrip_(files[4]);
               shadowfile = openac_(&ierr, shadowfilename, rwa);
               for (I = 0; I < search_.nstot; I++){
                    K = I+1;
                    if (SHADOW[I][0] != 0){
                         outputline = shdout_(&K, IMTYPE+I, SHADOW[I], 
                                 &tune4_.npar, which_model[I]);
                         fputs(outputline, shadowfile);
                    }
                    else{
                         outputline = shdout_(&K, IMTYPE+I, STARPAR[I],
                                 &tune4_.npar, which_model[I]);
                         fputs(outputline, shadowfile);
                    }
               }

               //also write errors on fits
               errfilename = stringstrip_(files[8]);
               errfile = openac_(&ierr, errfilename, rwa);
               for (I = 0; I < search_.nstot; I++){
                    K = I+1;
                    outputline = shdout_(&K, IMTYPE+I, SHADERR[I],
                                 &tune4_.npar, which_model[I]);
                    fputs(outputline, errfile);
               }

               fclose(shadowfile);
               free(shadowfilename);
               fclose(errfile);
               free(errfilename);
          }

          /* write covariance matrix file if flagged */
          if (flags[14][0] == 'Y'){
               covarfilename = stringstrip_(files[14]);
               covarfile = openac_(&ierr, covarfilename, rwa);
               for (I = 0; I < search_.nstot; I++){
                    K = I+1;
                    sprintf(outputline,"%4d %2d\n", K, IMTYPE[I]);
                    fputs(outputline, covarfile);
                    if (SHADOW[I][0] != 0){
                         for(J = 0; J < tune4_.nfit2; J++){
                              outputline = covarout_(SHADCOVAR[I][J]); 
                              fputs(outputline, covarfile);
                         }
                    }
               }
               fclose(covarfile);
               free(covarfilename);
          }

          /* write output image if flagged */
          if(flags[4][0] == 'Y'){ //output image
               newfits_(nx, ny, big, files[1], 1, files[0]);
          }
   
          if (first){
               first = 0;
          }

     }// end of while loop

     // make postage stamp (PS) images of all objects detected and residuals (RPS)
     // if flagged.  Stamps of the models (MPS) and images-neighbors (CPS)
     // are made in addstar called by improve in the while loop
     int len_pf_names = 0, len_rf_names = 0;
     char **pf_names, **rf_names ;
     char *pf_root, *rf_root, *dir;
     int indx_i;
     int ps_itype = 1;
     int** untouched_big;
     int dum_nx, dum_ny;
 
     if ((tune16_.flags[10][0] == 'Y') || (tune16_.flags[10][0] == 'y')){
          untouched_big = readfits_(files[0], &dum_nx, &dum_ny);          
          dir     = tune16_.files[9];
          pf_root = tune16_.files[10];
          len_pf_names = strlen(pf_root) + strlen(dir) + 10;
          pf_names = malloc_char_arr(search_.nstot, len_pf_names);
          for (indx_i = 0; indx_i < search_.nstot; indx_i++){
               sprintf(pf_names[indx_i],"%sd%04d%s.fits", dir, indx_i + 1, pf_root);
          }
          for (indx_i = 0; indx_i < search_.nstot; indx_i++){
               if ((IMTYPE[indx_i] == 1) && (tune22_.emenab)){
                    ps_itype = 10;
               }
               else{
                    if ((IMTYPE[indx_i] == 2) || (IMTYPE[indx_i] == 12)){
                         ps_itype = 2;
                    }
                    else{
                         ps_itype = 1;
                    }
               }
               stampout_(untouched_big, &nfast, &nslow, STARPAR[indx_i], 
                         ADDAREA[indx_i], ps_itype, pf_names[indx_i]);
          }
          free_char_arr(search_.nstot, pf_names);
          free_int_2darr(nslow, untouched_big);
     }
     if ((tune16_.flags[13][0] == 'Y') || (tune16_.flags[13][0] == 'y')){
          dir     = tune16_.files[9];
          rf_root = tune16_.files[13];
          len_rf_names = strlen(rf_root) + strlen(dir) + 10;
          rf_names = malloc_char_arr(search_.nstot, len_rf_names);
          for (indx_i = 0; indx_i < search_.nstot; indx_i++){
               sprintf(rf_names[indx_i],"%sd%04d%s.fits", dir, indx_i + 1, rf_root);
          }
          for (indx_i = 0; indx_i < search_.nstot; indx_i++){
               if ((IMTYPE[indx_i] == 1) && (tune22_.emenab)){
                    ps_itype = 10;
               }
               else{
                    if ((IMTYPE[indx_i] == 2) || (IMTYPE[indx_i] == 12)){
                         ps_itype = 2;
                    }
                    else{
                         ps_itype = 1;
                    }
               }
               stampout_(big, &nfast, &nslow, STARPAR[indx_i], 
                         ADDAREA[indx_i], ps_itype, rf_names[indx_i]);
          }
          free_char_arr(search_.nstot, rf_names);
     }

     free_int_2darr(nslow, big);
     free_int_2darr(nslow, noise);
       
     /* Close log file */
     if ((lverb > 1) && (strncmp(files[5], "TERM", 4) != 0)){
          closelog_();
     }
     endtime = clock();
     double runtime = (endtime-starttime)/(double)CLOCKS_PER_SEC;
     printf("Run time is %11g seconds\n", runtime);

     return 1;

}



