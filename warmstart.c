#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "logh.h"
#include "tuneable.h"
#include "search_struct.h"
#include "starlist_struct.h"
#include "addobj_struct.h"
#include "galpass_struct.h"
#include "cast_arr.h"
#include "io_dophotc.h"
#include "oblit.h"
#include "stdinpt.h"
#include "suminpt.h"
#include "addstar.h"
#include "warmstart.h"

/* dophot subroutine converted to c void functoin 03-30-2012 */

// WARNING not yet tested with a shadow file

//onestar is a pointer to the function name ONESTAR which returns double
void warmstart_( double (*ONESTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, char* filename )
{
     /* dereference pointers */
     int NFAST  = *NFAST_ptr;
     int NSLOW  = *NSLOW_ptr;

     /* rename used common block variables 
        which AREN'T changed frequently in subroutines */
     int NSTOT = search_.nstot;
     int*      IMTYPE = starlist_.IMTYPE;
     float**  STARPAR = starlist_.starpar;
     float**  SHADOW  = starlist_.shadow;
     float**  SHADERR = starlist_.shaderr;
     short int** ADDAREA = addobj_.addarea;

     int NPAR   = tune4_.npar;
     int lverb  = tune14_.lverb;
     float* AVA = tune15_.ava;
     char** files = tune16_.files;
     char** flags = tune16_.flags;

     /* substance of subroutine begins here */
     int I, J;
     static int first = 1;
     static float SUM[NPMAX];
     if (first){
          for(I = 0; I < NPMAX; I++){
               SUM[I] = 0;
          }
          first = 0;
     }
     char* prompt;
     FILE* objectfile;
     FILE* shadowfile;
     int mll = 300;
     char objectline[300];
     char shadowline[300];
     //object line set with sscanf in stdinpt, but can't hurt 
     char* rwa = "r"; //read toggle for open[sa]c
     int ISUB = -1;
     int DONE;
     int read_err;
     int IOSI, IOSC; //for reporting scanf and gets errors
     int NSTAR;
     int I1, I2, J1, J2;
     float SKY, TEMP;

     // open the object file.  if the given one is invalid, request a
     // new one
     objectfile = openac_(&read_err, files[2], rwa);
     if (read_err == 1){
          printf("Input object file doesn\'t exist on warmstart\n");
          prompt = "Enter input object file name";
          objectfile = opensc_(prompt, rwa);
     }

     // open the shadow file.  if the given one is invalid, request a
     // new one
     if (flags[6][0] == 'Y'){
          shadowfile = openac_(&read_err, files[6], rwa);
          if (read_err == 1){
               printf("Desired input shadow file doesn\'t exist on warmstart\n");
               prompt = "Enter input shadow file name";
               shadowfile = opensc_(prompt, rwa);
          }
     }
          
     IOSI  = 1; //all numbers read, no errors
     IOSC  = 1; //1 full line read, no errors
     I     = 0;
     NSTAR = 0;
     while ((IOSI == 1) && (IOSC == 1)){
//          printf("I, IOSI, IOSC = %d, %d, %d \n", I, IOSI, IOSC); 
          if (strncmp(flags[7], "INTER", 5) == 0){
               IOSC = (fgets(objectline, mll, objectfile) != NULL);
               IOSI = suminpt_(&I1, &I2, STARPAR[I], &NPAR, objectline);
          }
          if (strncmp(flags[7], "COMPL", 5) == 0){
               IOSC = (fgets(objectline, mll, objectfile) != NULL);
               IOSI = stdinpt_(&I1, &I2, STARPAR[I], &NPAR, objectline);
          }

          //continue only if lines read in properly from warmstart
          if ((IOSI == 1) && (IOSC == 1)){
               if (flags[6][0] == 'Y'){
                    IOSC = (fgets(shadowline, mll, shadowfile) != NULL);
                    IOSI = suminpt_(&J1, &J2, SHADOW[I], &NPAR, shadowline);
                    SKY = SHADOW[I][0];

                    // Changed indices.
                    if (SHADOW[I][1] > 0.0f){
                         TEMP = (SHADOW[I][1] + SKY)*
                                  (1.0f/(SHADOW[I][1]*SHADOW[I][1]));
                    }
                    else{
                         TEMP = 1.0E10f;
                    }
                    for(J = 0; J < NPAR; J++){
                         SHADERR[I][J] = TEMP;
                    }

                    if ((J1 != I1) || (J2 != I2)){
                         fprintf(logfile,"Trouble on WARMSTART:\n");
                         fprintf(logfile,"Object_input and Shadow_input ");
                         fprintf(logfile,"files have discrepant data!\n");
                         fprintf(logfile,"Forcing STOP!\n");
                         fprintf(logfile,"J1 != I1 OR J2 != I2 , %d %d %d %d\n",
                                          J1, I1, J2, I2);
                         exit(6666);
                    } 
               }
 
               NSTOT = I+1;
               IMTYPE[I] = I2;
               if ((IMTYPE[I] == 1) || (IMTYPE[I] == 11)){
                    NSTAR += 1;
                    for(J = 0; J < NPAR; J++){
                         SUM[J] += STARPAR[I][J];
                    } 
               }
               else{
                    if (IMTYPE[I] == 8){
                         if ((STARPAR[I][5] == 0.0f ) ||
                             (STARPAR[I][5] == -1.0f) || 
                             (STARPAR[I][5] == 90.0f)) {
                              DONE = oblit_(ONESTAR, BIG, NOISE,
                                            &NFAST, &NSLOW, STARPAR[I]);
                         }
                         else{
                              DONE = oblit_ellipse_(ONESTAR, BIG, NOISE,
                                            &NFAST, &NSLOW, STARPAR[I]);
                         }
                    }
               }

               if ((IMTYPE[I] != 6) && (IMTYPE[I] != 8)){
                    if (IMTYPE[I] == 2) {
                         galpass_.bigfoot = 1;
                    }
                    addstar_(ONESTAR, BIG, NOISE, NFAST, NSLOW,
                                      STARPAR[I], ADDAREA[I], ISUB, 
                                      0, " ", 0, " ");
                    if (IMTYPE[I] == 2) {
                         galpass_.bigfoot = 0;
                    }
               }
          } //end IOS good if
          I += 1;
     } //end IOS good while

     if (lverb > 10){
          fprintf(logfile,"\n");
          fprintf(logfile,"%d objects subtracted from image at WARMSTART\n",
                           NSTOT);
     }

     if (NSTAR != 0){
          for(I = 0; I < NPAR; I++){
               AVA[I] = SUM[I]/(float)(NSTAR);
          }
     }

     /* closing files */
     fclose(objectfile);
     if (flags[6][0] == 'Y'){
          fclose(shadowfile);
     }

     /* recasting changed common blocks */
     search_.nstot = NSTOT;

}
