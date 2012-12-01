#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tuneable.h"
#include "mini_mathlib.h"
#include "cast_arr.h"
#include "io_dophotc.h"
#include "paramfile.h"
#include "tuneup.h"

/* dophot subroutine converted to  c void function 04-01-2012 */

void tuneup_(int command_line_pm, char* pm_file_name)
{
     int I; //iterator here forward

//printf("here0\n");
     // get the default and modified parameters from the files
     char** keywords = malloc(nphead*sizeof(char*));
     char** items    = malloc(nphead*sizeof(char*));
     char** comments = malloc(nphead*sizeof(char*));
     int nlines, nlines_init;

//printf("here1\n");
     paramfile_(command_line_pm, pm_file_name, keywords, items, comments, &nlines);
     nlines_init = nlines; // used for keeping track of how many keywords
                           // and comments ultimately added

//printf("here2\n");
     /* Read the file names and file usage flags.  These variables
        will be transmitted to the outside world in the array FILES. */
     char** files = malloc(NFF*sizeof(char*)); // individual str malloced
     char** flags = malloc(NFF*sizeof(char*)); // individually to proper len

//printf("here3\n");
     /* Read the input image name.  THE NAME MUST EXIST.  If not, it is
        requested. */
     FILE* tempfile;
     files[0] = get_string_item_(keywords, items, "IMAGE_IN", nlines);
     /* check if file will open. if file wont open, request a valid file */
     int readerr = 0;
     char* rwa = "r";
     char* prompt;
     char tempstring[80];
     int still_needs_file;
     tempfile = openac_(&readerr, files[0], rwa);
     if (readerr == 1){
          free(files[0]); //in case length is very small
          printf("Input image does not exist\n");
          prompt = "Enter input image name";
          //NOTE:  I'm essentially rewritting opensc here because I want to
          //       keep the image name in files[0]
          still_needs_file = 1; //true
          while (still_needs_file){
               printf("%s : ", prompt);
               if (scanf("%s", tempstring) != 1){
                    printf("Error in file name.  Try again.\n");
               }
               else{
                    tempfile = fopen(tempstring, rwa);
                    if (tempfile == NULL){
                         printf("Error opening file.  Try again.\n");
                    }
                    else{
                         still_needs_file = 0; //false, exit routine
                    }
               }
          }
          files[0] = malloc((strlen(tempstring)+1)*sizeof(char));
          strcpy(files[0], tempstring);
          put_item_(keywords, items, "IMAGE_IN", &nlines, files[0], 2);
     }
     fclose(tempfile);
          
//printf("here4\n");
     /* Determine if an output picture file is to be saved. */
     flags[4] = malloc(4*sizeof(char));
     files[1] = get_string_item_(keywords, items, "IMAGE_OUT", nlines);
     if(files[1][0] == '\0'){ //no file
          strcpy(flags[4], "NO ");
     }
     else{
          strcpy(flags[4], "YES");
     }
     
//printf("here5\n");
     /* Read input objects list.  Leave null if none is given. */
     flags[5] = malloc(4*sizeof(char));
     files[2] = get_string_item_(keywords, items, "OBJECTS_IN", nlines);
     if(files[2][0] == '\0'){ //no file
          strcpy(flags[5], "NO ");
     }
     else{
          strcpy(flags[5], "YES");
     }

//printf("here6\n");
     /* The output objects file name MUST EXIST.  If it is not in the 
        parameter file, it is requested. */
     files[3] = get_string_item_(keywords, items, "OBJECTS_OUT", nlines);
     if(files[3][0] == '\0'){ //no file
          free(files[3]);
          printf("Enter output object file name: ");
          scanf("%s",tempstring);
          files[3] = malloc((strlen(tempstring)+1)*sizeof(char));
          put_item_(keywords, items, "OBJECTS_OUT", &nlines, files[3], 2);
     }

//printf("here7\n");
     /* Ask for the shadow file name if it is to be saved.  If no name 
        is provided, the output file will get a default designation. */
     flags[3] = malloc(4*sizeof(char));
     files[4] = get_string_item_(keywords, items, "SHADOWFILE_OUT", nlines);
     if(files[4][0] == '\0'){ //no file
          strcpy(flags[3], "NO ");
     }
     else{
          strcpy(flags[3], "YES");
     }

//printf("here8\n");
     /* Figure out about an input shadow file if desired. */
     flags[6] = malloc(4*sizeof(char));
     if (strcmp(flags[5], "YES") == 0){ //if there are objects_in
          files[6] = get_string_item_(keywords, items, "SHADOWFILE_IN", nlines);
          if(files[6][0] == '\0'){ //no file
               strcpy(flags[6], "NO ");
          }
          else{
               strcpy(flags[6], "YES");
          }
     }
     else{
          files[6] = malloc(1*sizeof(char));
          files[6][0] = '\0';
          strcpy(flags[6], "NO ");
     }

//printf("here9\n");
     /* Ask for the log file name if it is to be saved.  If no name 
        is provided, the output file will get a default designation. */
     tune14_.lverb = get_int_item_(keywords, items, "LOGVERBOSITY", nlines);
     if (tune14_.lverb > 0){
          files[5] = get_string_item_(keywords, items, "LOGFILE", nlines);
          if(files[5][0] == '\0'){ //no file
               free(files[5]);
               files[5] = malloc((strlen("logfile.dat")+1)*sizeof(char));
               strcpy(files[5],"logfile.dat");
               put_item_(keywords, items, "LOGFILE", &nlines, files[5], 2);
          }
     }
     else{
          files[5] = malloc(1*sizeof(char));
          files[5][0] = '\0';
     }
     
//printf("here10\n");
     /* Now read the psf type flag.  If no parameters are present, the
        defaults PGAUSS is assumed.  infer the number of fit parameters
        from the model type */
     flags[0] = get_string_item_(keywords, items, "PSFTYPE", nlines);
     if (flags[0][0] == '\0'){ //no psf type given
          free(flags[0]);
          flags[0] = malloc((strlen("PGAUSS")+1)*sizeof(char));
          strcpy(flags[0],"PGAUSS");
          put_item_(keywords, items, "PSFTYPE", &nlines, flags[0], 2);
     }
     else{
          if ( (strncmp(flags[0], "PGAUSS", 5) != 0) &&
               (strncmp(flags[0], "GAUSS", 5)  != 0) && 
               (strncmp(flags[0], "PGALTMODEL", 5) != 0) && 
               (strncmp(flags[0], "GALTMODEL", 5)  != 0) && 
               (strncmp(flags[0], "SERSIC", 5) != 0) && 
               (strncmp(flags[0], "EXTPGAUSS", 5) != 0) ){

               free(flags[0]);
               do{
                    printf("Enter valid model type (PGAUSS, GAUSS): "); 
                    scanf("%s",tempstring);
               }while( (strncmp(tempstring, "PGAUSS", 5) != 0) &&
                       (strncmp(tempstring, "GAUSS", 5)  != 0) &&
                       (strncmp(tempstring, "PGALTMODEL", 5) != 0) &&
                       (strncmp(tempstring, "GALTMODEL", 5)  != 0) &&
                       (strncmp(tempstring, "SERSIC", 5) != 0) && 
                       (strncmp(tempstring, "EXTPGAUSS", 5) != 0) );
          
               flags[0] = malloc((strlen(tempstring)+1)*sizeof(char));
               strcpy(flags[0], tempstring);
               put_item_(keywords, items, "PSFTYPE", &nlines, flags[0], 2);
          }
     }
     tune4_.npar  = 7; //NPARAM
     tune4_.nfit2 = 7; //NFITSHAPE
     if((strncmp(flags[0], "GALTMODEL",  5) == 0) ||
        (strncmp(flags[0], "PGALTMODEL", 5) == 0)){
          tune4_.npar  += 4;
          tune4_.nfit2 += 4;
     }
     if(strncmp(flags[0], "SERSIC", 5) == 0){
          tune4_.npar  += 1;
          tune4_.nfit2 += 1;
     }
     if(strncmp(flags[0], "EXTPGAUSS", 5) == 0){
          tune4_.npar  += 2;
          tune4_.nfit2 += 2;
     }

//printf("here11\n");
     /* Now read the sky type flag.  If no parameters are present, the
        defaults PLANE is assumed.  If an invalid type is given, a
        valid type is requested */
     flags[1] = get_string_item_(keywords, items, "SKYTYPE", nlines);
     if (flags[1][0] == '\0'){ //no sky type given
          free(flags[1]);
          flags[1] = malloc((strlen("PLANE")+1)*sizeof(char));
          strcpy(flags[1],"PLANE");
          put_item_(keywords, items, "SKYTYPE", &nlines, flags[1], 2);
     }
     if ( (strncmp(flags[1], "PLANE", 5) != 0) &&
          (strncmp(flags[1], "MEDIA", 5) != 0) &&
          (strncmp(flags[1], "HUBBL", 5) != 0) ){
          free(flags[1]);
          do{
               printf("Enter valid sky type (PLANE, MEDIAN, HUBBLE): "); 
               scanf("%s",tempstring);
          }while( (strncmp(tempstring, "PLANE", 5) != 0) &&
                  (strncmp(tempstring, "MEDIA", 5) != 0) &&
                  (strncmp(tempstring, "HUBBL", 5) != 0) );
          flags[1] = malloc((strlen(tempstring)+1)*sizeof(char));
          strcpy(flags[1], tempstring);
          put_item_(keywords, items, "SKYTYPE", &nlines, flags[1], 2);
     }

//printf("here12\n");
     /* get type of input object file.  if none given, assume INTERNAL */
     /* if invalid type, request valid type */
     flags[7] = get_string_item_(keywords, items, "OBJTYPE_IN", nlines);
     if (flags[7][0] == '\0'){ //object type given
          free(flags[7]);
          flags[7] = malloc((strlen("INTERNAL")+1)*sizeof(char));
          strcpy(flags[7],"INTERNAL");
          put_item_(keywords, items, "OBJTYPE_IN", &nlines, flags[7], 2);
     }
     if ( (strncmp(flags[7], "COMPL", 5) != 0) &&
          (strncmp(flags[7], "INTER", 5) != 0) ){
          free(flags[7]);
          do{
               printf("Enter valid input object format type (COMPLETE, INTERNAL): "); 
               scanf("%s",tempstring);
          }while( (strncmp(tempstring, "COMPL", 5) != 0) &&
                  (strncmp(tempstring, "INTER", 5) != 0) );
          flags[7] = malloc((strlen(tempstring)+1)*sizeof(char));
          strcpy(flags[7], tempstring);
          put_item_(keywords, items, "OBJTYPE_IN", &nlines, flags[7], 2);
     }
     
//printf("here13\n");
     /* get type of output object file.  if none given, assume INTERNAL */
     /* if invalid type, request valid type */
     flags[2] = get_string_item_(keywords, items, "OBJTYPE_OUT", nlines);
     if (flags[2][0] == '\0'){ //object type given
          free(flags[2]);
          flags[2] = malloc(9*sizeof(char));
          strcpy(flags[2],"INTERNAL");
          put_item_(keywords, items, "OBJTYPE_OUT", &nlines, flags[2], 2);
     }
     if ( (strncmp(flags[2], "COMPL", 5) != 0) &&
          (strncmp(flags[2], "INCOM", 5) != 0) &&
          (strncmp(flags[2], "INTER", 5) != 0) ){
          free(flags[2]);
          do{
               strcpy(flags[2], "        ");
               printf("Enter valid output object format type \n");
               printf("(COMPLETE, INCOMPLETE, INTERNAL): "); 
               scanf("%s", tempstring);
          }while( (strncmp(tempstring, "COMPL", 5) != 0) &&
                  (strncmp(tempstring, "INCOM", 5) != 0) &&
                  (strncmp(tempstring, "INTER", 5) != 0) );
          flags[2] = malloc((strlen(tempstring)+1)*sizeof(char));
          strcpy(flags[2], tempstring);
          put_item_(keywords, items, "OBJTYPE_OUT", &nlines, flags[2], 2);
     }

//printf("files and flags read\n");
     /* get non-flag/file tuneable keywords from param files */

     /* tune1_. gxwid, gywid, skyguess*/
     float tilt       = get_float_item_(keywords, items, "TILT",       nlines);
     if( tilt == -999.9f) tilt = 0.0f;

     float axis_ratio = get_float_item_(keywords, items, "AXIS_RATIO", nlines);
     if( axis_ratio == -999.9f) axis_ratio = 0.0f;

     float fwhm       = get_float_item_(keywords, items, "FWHM",       nlines);
     if( fwhm == -999.9f){
         fwhm = 4.0f;
         printf("Problem with FWHM parameter.  Setting to 4.0 \n");
         printf("Update this parameter for better results.\n");
     }
     fwhm = fwhm*1.2f;

     tilt = tilt/57.29578f;
     float gmajwid = (fwhm/2.3548f)*(fwhm/2.3548f);
     tune1_.gxwid = gmajwid*
            (cosf(tilt)*cosf(tilt) + axis_ratio*sinf(tilt)*sinf(tilt));
     tune1_.gywid = gmajwid*
            (axis_ratio*cosf(tilt)*cosf(tilt) + sinf(tilt)*sinf(tilt));
     tune1_.gxwid = (float)( (int)(10.0f*tune1_.gxwid + 0.5f) )/ 10.0f;
     tune1_.gywid = (float)( (int)(10.0f*tune1_.gywid + 0.5f) )/ 10.0f;
     float fwhmx = 2.3548f*sqrtf(tune1_.gxwid);
     float fwhmy = 2.3548f*sqrtf(tune1_.gywid);
     if (tune14_.lverb > 10){
          printf("fwhmx, fwhmy = %f %f \n", fwhmx, fwhmy);
     }

     tune1_.skyguess  = get_float_item_(keywords, items, "SKY",       nlines);
     if( tune1_.skyguess == -999.9f){
         tune1_.skyguess = 0.0f;
         printf("Problem with SKY parameter.  Setting to 0.0 \n");
         printf("Update this parameter for better results.\n");
     }

//printf("here14\n");
     /* tune15_.ava */
     float* ava = malloc_float_1darr(NPMAX);
     ava[0] = tune1_.skyguess;
     ava[1] = 0.0f;
     ava[2] = 0.0f;
     ava[3] = 0.0f;
     ava[4] = tune1_.gxwid;
     ava[5] = 0.01f/sqrtf(tune1_.gywid*tune1_.gxwid);
     ava[6] = tune1_.gywid;
     if((strncmp(flags[0], "GALTMODEL",  5 ) == 0) ||
        (strncmp(flags[0], "PGALTMODEL", 5 ) == 0)){
          ava[7]  = 0.0f;
          ava[8]  = 0.0f;
          ava[9]  = 0.0f;
          ava[10] = 0.0f;
     }
     if(strncmp(flags[0], "SERSIC", 5) == 0) ava[7] = 2.0f; //sersic index
     if(strncmp(flags[0], "EXTPGAUSS", 5) == 0) {
//          ava[7] = 2.5f; //B4
          ava[7] = 1.0f; //B4
          ava[8] = 1.0f; //B6
     }
     tune15_.ava = ava;
     
//printf("here15\n");
     /* tune15_.ava */
     /* tune11_.eperdn, rdnoise */
     tune11_.eperdn = get_float_item_(keywords, items, "EPERDN", nlines);
     if( tune11_.eperdn == -999.9f) tune11_.eperdn = 1.0f;

     tune11_.rnoise = get_float_item_(keywords, items, "RDNOISE", nlines);
     if( tune11_.rnoise == -999.9f) tune11_.rnoise = 1.0f;

//printf("here16\n");
     /* Now determine if autoscaling is desired. */
     char* autoscale = get_string_item_(keywords, items, "AUTOSCALE", nlines);
     int do_autoscale;
     if ((autoscale[0] == '\0') ||
         (autoscale[0] == 'Y' ) || (autoscale[0] == 'y')){
          do_autoscale = 1; //default or YES autoscale
     }
     else{
          do_autoscale = 0; //if keyword is specified but not yes, then assumed no
     }
     free(autoscale);

//printf("here17\n");
     /* autoscale parameters */
     float scalefb, fbmin, scaleab, abmin, scalemb, ambmin;
     short int i4;
     short int* irect = malloc_si_1darr(2);
     float* arect     = malloc_float_1darr(2);
     char i4_char[81];
     char arect_char[81];
     char irby2_char[81];
     if (do_autoscale){
          scalefb = get_float_item_(keywords, items, "SCALEFITBOX",  nlines);
          fbmin   = get_float_item_(keywords, items, "FITBOXMIN",    nlines);
          scaleab = get_float_item_(keywords, items, "SCALEAPBOX",   nlines);
          abmin   = get_float_item_(keywords, items, "APBOXMIN",     nlines);
          scalemb = get_float_item_(keywords, items, "SCALEMASKBOX", nlines);
          ambmin  = get_float_item_(keywords, items, "AMASKBOXMIN",  nlines);
          //assumes these items will exist and not default to -999.9 

          /* Apply these factors to the relevant quantities. */
          /* tune2_.irect */
          i4 = (short int)(max1(fwhmx*scalefb,fbmin) + 0.5f);
          if ( (i4 % 2) == 0 ) i4 += 1;
          irect[0] = i4;
          sprintf(i4_char, "%d", i4);
          put_item_(keywords, items, "NFITBOX_X", &nlines, i4_char, 1);
          
          i4 = (short int)(max1(fwhmy*scalefb,fbmin) + 0.5f);
          if ( (i4 % 2) == 0 ) i4 += 1;
          irect[1] = i4;
          sprintf(i4_char, "%d", i4);
          put_item_(keywords, items, "NFITBOX_Y", &nlines, i4_char, 1);
          
          tune2_.irect = irect;

          /* tune2_.arect */
          arect[0] = (float)(int)(max1(fwhmx*scaleab,abmin) + 0.5f);
          if ( ((int)(arect[0]) % 2) == 0 ) arect[0] += 1.0f;
          sprintf(arect_char, "%f", arect[0]);
          put_item_(keywords, items, "APBOX_X", &nlines, arect_char, 0);
          
          arect[1] = (float)(int)(max1(fwhmy*scaleab,abmin) + 0.5f);
          if ( ((int)(arect[1]) % 2) == 0 ) arect[1] += 1.0f;
          sprintf(arect_char, "%f", arect[1]);
          put_item_(keywords, items, "APBOX_Y", &nlines, arect_char, 0);
          
          tune2_.arect = arect;

          /* tune12_.ixby2, iyby2 */          
          tune12_.ixby2 = (int)(max1(fwhmx*scalemb,ambmin) + 0.5f);
          if ( (tune12_.ixby2 % 2) == 0 ) tune12_.ixby2 += 1;
          tune12_.ixby2 = (tune12_.ixby2 - 1)/2;
          sprintf(irby2_char, "%d", tune12_.ixby2);
          put_item_(keywords, items, "MASKBOX_X", &nlines, irby2_char, 1);

          tune12_.iyby2 = (int)(max1(fwhmy*scalemb,ambmin) + 0.5f);
          if ( (tune12_.iyby2 % 2) == 0 ) tune12_.iyby2 += 1;
          tune12_.iyby2 = (tune12_.iyby2 - 1)/2;
          sprintf(irby2_char, "%d", tune12_.iyby2);
          put_item_(keywords, items, "MASKBOX_Y", &nlines, irby2_char, 1);

     }
     else{ //if not autoscale, read parameters directly from file     
          /* tune2_.irect, arect */
          irect[0] = (short int)get_int_item_(keywords, items, "NFITBOX_X", nlines);
          irect[1] = (short int)get_int_item_(keywords, items, "NFITBOX_Y", nlines);
          tune2_.irect = irect;

          arect[0] = get_float_item_(keywords, items, "APBOX_X",   nlines);
          arect[1] = get_float_item_(keywords, items, "APBOX_Y",   nlines);
          tune2_.arect = arect;
 
          /* tune12_.ixby2, iyby2 */          
          tune12_.ixby2 = get_int_item_(keywords, items, "MASKBOX_X", nlines);
          if ( (tune12_.ixby2 % 2) == 0 ) tune12_.ixby2 += 1;
          tune12_.ixby2 = (tune12_.ixby2 - 1)/2;

          tune12_.iyby2 = get_int_item_(keywords, items, "MASKBOX_Y", nlines);
          if ( (tune12_.iyby2 % 2) == 0 ) tune12_.iyby2 += 1;
          tune12_.iyby2 = (tune12_.iyby2 - 1)/2;

//          printf("irect = (%d, %d)\n", tune2_.irect[0], tune2_.irect[1]); 
//          printf("arect = (%f, %f)\n", tune2_.arect[0], tune2_.arect[1]); 
//          printf("i by2 = (%d, %d)\n", tune12_.ixby2, tune12_.iyby2); 
     }

//printf("here18\n");
     /* Ask if positions will be fixed. tune21_.fixpos (int)*/
     char* fixpos_yn = get_string_item_(keywords, items, "FIXPOS", nlines);
     if ((fixpos_yn[0] == 'Y' ) || (fixpos_yn[0] == 'y')){
          tune21_.fixpos = 1; //YES fix position
     }
     else{
          tune21_.fixpos = 0; //if keyword is not specified or not yes, then no
     }
     free(fixpos_yn);

     
//printf("here19\n");
     /* Ask if positions will be fixed. tune21_.fixpos (int)*/
     int do_autothresh;
     char* autothresh = get_string_item_(keywords, items, "AUTOTHRESH", nlines);
     if ((autothresh[0] == 'Y' ) || (autothresh[0] == 'y')){
          do_autothresh = 1; //YES fix position
     }
     else{
          do_autothresh = 0; //if keyword is not specified or not yes, then no
     }
     free(autothresh);

//printf("here20\n");
     /* read in autothresh parameters if do_autothresh */
     float sigbot, sigthresh, sigskyguess;
     if (do_autothresh){
          sigbot    = get_float_item_(keywords, items, "SIGMAIBOTTOM",   nlines);
          sigthresh = get_float_item_(keywords, items, "SIGMATHRESHMIN", nlines);
          sigskyguess = sqrtf(  tune11_.eperdn*tune1_.skyguess 
                              + tune11_.rnoise*tune11_.rnoise  )/tune11_.eperdn;
          
          /* Now, calculate the sigma of the sky and apply these parameters. */
          /* tune3_. ibot, tmin */ 
          tune3_.ibot = (int)(tune1_.skyguess - sigbot*sigskyguess);
          sprintf(tempstring, "%d", tune3_.ibot);
          put_item_(keywords, items, "IBOTTOM", &nlines, tempstring, 1);

          tune3_.tmin = sigthresh*sigskyguess;
          sprintf(tempstring, "%f", tune3_.tmin);
          put_item_(keywords, items, "THRESHMIN", &nlines, tempstring, 0);
     }
     else{ // if not auto, read parameters in directly from file
          tune3_.ibot =   get_int_item_(keywords, items, "IBOTTOM",   nlines);
          tune3_.tmin = get_float_item_(keywords, items, "THRESHMIN", nlines);
     }
          
//printf("here21\n");
     /* tune18_.nthpix (int)*/
     tune18_.nthpix = 0; 
     /* if nthpix is carried as 0 DoPHOT will change it to sqrt(nfast)
        nfast is not known till data file is read, so it cannot be set 
        here. */

//printf("here22\n");
     /* set MEDIAN sky keywords, tune20_, tune18_.nthpix */
     if (strncmp(flags[1], "MEDIA", 5) == 0){
          tune20_.jhxwid = get_int_item_(keywords, items, "JHXWID", nlines);
          tune20_.jhywid = get_int_item_(keywords, items, "JHYWID", nlines);
          tune20_.mprec  = get_int_item_(keywords, items, "MPREC",  nlines);
          tune18_.nthpix = get_int_item_(keywords, items, "NTHPIX",  nlines);
          if (tune20_.jhxwid <= 0){
               tune20_.jhxwid = max(2,(int)(4.0f*sqrtf(tune1_.gxwid) + 0.5f));
          }
          if (tune20_.jhywid <= 0){
               tune20_.jhywid = max(2,(int)(4.0f*sqrtf(tune1_.gywid) + 0.5f));
          }
          if (tune20_.mprec  <= 0){ 
               tune20_.mprec  = (int)(tune3_.tmin/4.0f);
          }
          tune20_.mprec = max(1, tune20_.mprec);
     }

//printf("here23\n");
     /* Now, read all remaining numeric parameters.  There's still a bunch. */
     short int* krect = malloc_si_1darr(2);
     krect[0] = (short int)get_int_item_(keywords, items, "NFITBOXFIRST_X", nlines);
     krect[1] = (short int)get_int_item_(keywords, items, "NFITBOXFIRST_Y", nlines);
     tune2_.krect = krect;

     tune3_.itop     =   get_int_item_(keywords, items, "ITOP", nlines);
     tune3_.tmax     = get_float_item_(keywords, items, "THRESHMAX", nlines);
     tune3_.tfac     = get_float_item_(keywords, items, "THRESHDEC", nlines);

     tune4_.nit      =   get_int_item_(keywords, items, "NFITITER", nlines);
//     tune4_.npar     =   get_int_item_(keywords, items, "NPARAM", nlines);
     tune4_.nfit1    =   get_int_item_(keywords, items, "NFITMAG", nlines);
//     tune4_.nfit2    =   get_int_item_(keywords, items, "NFITSHAPE", nlines);

     tune5_.gxpnd    = 3.0f;// (float)  KLUGE
     tune5_.fac      = get_float_item_(keywords, items, "RESIDNOISE", nlines);
     tune5_.xpnd     = get_float_item_(keywords, items, "FOOTPRINT_NOISE", nlines);

     tune6_.nphsub   =   get_int_item_(keywords, items, "NPHSUB", nlines);
     tune6_.nphob    =   get_int_item_(keywords, items, "NPHOB", nlines);
     tune6_.icrit    =   get_int_item_(keywords, items, "ICRIT", nlines);
     tune6_.cmax     = get_float_item_(keywords, items, "CENTINTMAX", nlines);
     tune6_.ctpersat = get_float_item_(keywords, items, "CTPERSAT", nlines);
     tune6_.widobl   = get_float_item_(keywords, items, "COSOBLSIZE", nlines);

     tune7_.n0left   = 0  ;// (int) KLUGE
     tune7_.n0right  = 0  ;// (int) KLUGE

     tune8_.stograt  = get_float_item_(keywords, items, "STARGALKNOB", nlines);
     tune8_.discrim  = get_float_item_(keywords, items, "STARCOSKNOB", nlines);
     float* sig = malloc_float_1darr(3);
     sig[0]          = get_float_item_(keywords, items, "SIGMA1", nlines);
     sig[1]          = get_float_item_(keywords, items, "SIGMA2", nlines);
     sig[2]          = get_float_item_(keywords, items, "SIGMA3", nlines);
     tune8_.sig = sig;

     tune9_.crit7    = get_float_item_(keywords, items, "SNLIM7", nlines);
     tune9_.snlim    = get_float_item_(keywords, items, "SNLIM", nlines);
     tune9_.bumpcrit = get_float_item_(keywords, items, "SNLIMMASK", nlines);
     tune9_.sn2cos   = get_float_item_(keywords, items, "SNLIMCOS", nlines);
     tune9_.chicrit  = get_float_item_(keywords, items, "CHI2MINBIG", nlines);
     tune9_.xtra     = get_float_item_(keywords, items, "XTRA", nlines);
     tune9_.crit7    = tune9_.crit7*tune9_.crit7;
     tune9_.sn2cos   = tune9_.sn2cos*tune9_.sn2cos;

     tune10_.enuff4  = get_float_item_(keywords, items, "ENUFF4", nlines);
     tune10_.enuff7  = get_float_item_(keywords, items, "ENUFF7", nlines);

     tune13_.nbadleft  = get_int_item_(keywords, items, "NBADLEFT", nlines);
     tune13_.nbadright = get_int_item_(keywords, items, "NBADRIGHT", nlines);
     tune13_.nbadtop   = get_int_item_(keywords, items, "NBADTOP", nlines);
     tune13_.nbadbot   = get_int_item_(keywords, items, "NBADBOT", nlines);

     tune15_.max_perf = get_int_item_(keywords, items, "MAX_PERF", nlines);

     float* acc = malloc_float_1darr(NPMAX);
     acc[0]          = get_float_item_(keywords, items, "RELACC1", nlines);
     acc[1]          = get_float_item_(keywords, items, "RELACC4", nlines);
     acc[2]          = get_float_item_(keywords, items, "RELACC2", nlines);
     acc[3]          = get_float_item_(keywords, items, "RELACC3", nlines);
    
     char relacc_str[80]; 
     for (I = 4; I < tune4_.nfit2; I++){
          sprintf(relacc_str, "RELACC%d", I+1);
          acc[I] = get_float_item_(keywords, items, relacc_str, nlines);
     } 
     // NOTE INTENTIONAL WRONG ORDER OF RELACC & ALIM FOR HISTORICAL REASONS
     tune15_.acc = acc;

     float* alim = malloc_float_1darr(NPMAX);
     alim[0]         = get_float_item_(keywords, items, "ABSLIM1", nlines);
     alim[1]         = get_float_item_(keywords, items, "ABSLIM4", nlines);
     alim[2]         = get_float_item_(keywords, items, "ABSLIM2", nlines);
     alim[3]         = get_float_item_(keywords, items, "ABSLIM3", nlines);
     alim[4]         = get_float_item_(keywords, items, "ABSLIM5", nlines);
     alim[5]         = get_float_item_(keywords, items, "ABSLIM6", nlines);
     alim[6]         = get_float_item_(keywords, items, "ABSLIM7", nlines);
     if(strncmp(flags[0], "SERSIC", 5) == 0){
          alim[7] = -10.0f;
     }
     if(strncmp(flags[0], "EXTPGAUSS", 5) == 0){
          alim[7] = -10.0f;
          alim[8] = -10.0f;
     }
     if((strncmp(flags[0], "GALTMODEL",  5 ) == 0) ||
        (strncmp(flags[0], "PGALTMODEL", 5 ) == 0)){
          alim[7]  = -1.0f;
          alim[8]  = -1.0f;
          alim[9]  = -1.0f;
          alim[10] = -1.0f;
     }
     tune15_.alim = alim;

     tune17_.beta4   = get_float_item_(keywords, items, "BETA4", nlines);
     tune17_.beta6   = get_float_item_(keywords, items, "BETA6", nlines);

     tune18_.pixthresh   = get_float_item_(keywords, items, "PIXTHRESH", nlines);

     tune19_.apmagmaxerr = get_float_item_(keywords, items, "APMAG_MAXERR", nlines);


     /* BEGIN MODIFICATION FOR SCHECHTER'S EMPIRICAL FITTING */     

//printf("here24\n");
     /* Determine if empirical  subraster is to be saved. */
     flags[8] = malloc(4*sizeof(char));
     files[7] = get_string_item_(keywords, items, "EMP_SUBRAS_OUT", nlines);
     if (files[7][0] == '\0'){// no file specified
          strcpy(flags[8], "NO ");
     }
     else{ 
          strcpy(flags[8], "YES");
     }
     
     tune22_.ihtab   = get_int_item_(keywords, items, "NEMP_PSF_BOX", nlines);
     tune22_.ihtab  -= 1;
     tune22_.ihtab  -= tune22_.ihtab/2;
     tune22_.ihtab   = min(tune22_.ihtab, IHSIDE);
     /* If even forces odd box size, one pixel smaller. */

     tune22_.eminsig = get_float_item_(keywords, items, "EMP_REJ_RAD_SIG", nlines);
     tune22_.efac    = get_float_item_(keywords, items, "EMP_RESIDNOISE", nlines);
     tune22_.emthrsh = get_float_item_(keywords, items, "THRESHEMP", nlines);
     tune22_.gfac    = 1.0f;// (float) KLUGE

//printf("here25\n");
     /* X and Y positions and Z value for empirical subraster */
     tune23_.nempski = get_int_item_(keywords, items, "N_EMP_SKIP", nlines);
     flags[9] = malloc(4*sizeof(char));
     tune23_.targx = get_float_item_(keywords, items, "EMP_STAR_X", nlines);
     tune23_.targy = get_float_item_(keywords, items, "EMP_STAR_Y", nlines);
     tune23_.targz = get_float_item_(keywords, items, "EMP_STAR_Z", nlines);
     int chooemp;
     chooemp = ((tune23_.targx > 0) && (tune23_.targy > 0) && (tune23_.targz > 0));
     if (chooemp){ //empirical star specified x,y,and z
          strcpy(flags[9], "YES");
     }
     else{ 
          strcpy(flags[9], "NO ");
     }

     /* END SCHECHTER'S MODIFICATIIONS */

     /* BEGIN LEVINSON MODIFICATION FOR ERRORS_OUT &
        COVARIANCE MATRICES files */
     /* if shadow files requested, also output error file */
     if (strncmp(flags[3], "YES", 3) == 0){
          files[8] = get_string_item_(keywords, items, "ERRORS_OUT", nlines);
     } else{
          files[8] = malloc(1*sizeof(char));
          files[8][0] = '\0';
     }

     /* Covariance matrix file */
     flags[14] = malloc(4*sizeof(char));
     files[14] = get_string_item_(keywords, items, "COVARS_OUT", nlines);
     if(files[14][0] == '\0'){ //no file
          strcpy(flags[14], "NO ");
     }
     else{
          strcpy(flags[14], "YES");
     }
     /* END LEVINSON MODIFICATION FOR ERRORS_OUT & COVARS_OUT files */

     /* BEGIN LEVINSON MODIFICATION FOR POSTAGE STAMP MODEL FILE OUTPUTS */ 

     /* Determine if postage stamp files for each object are to be saved */
     // so flags[10-13] align with files[10-13]
     // POSTAGE STAMP (PS) from original image
     flags[10] = malloc(4*sizeof(char));
     files[10] = get_string_item_(keywords, items, "PSNAME_ROOT", nlines);
     if(files[10][0] == '\0'){ //no file
          strcpy(flags[10], "NO ");
     }
     else{
          strcpy(flags[10], "YES");
     }

     // PS from Cleaned, neighbor subtracted image 
     flags[11] = malloc(4*sizeof(char));
     files[11] = get_string_item_(keywords, items, "CPSNAME_ROOT", nlines);
     if(files[11][0] == '\0'){ //no file
          strcpy(flags[11], "NO ");
     }
     else{
          strcpy(flags[11], "YES");
     }
    
     // PS of Model image 
     flags[12] = malloc(4*sizeof(char));
     files[12] = get_string_item_(keywords, items, "MPSNAME_ROOT", nlines);
     if(files[12][0] == '\0'){ //no file
          strcpy(flags[12], "NO ");
     }
     else{
          strcpy(flags[12], "YES");
     }
    
     // PS of Residual image 
     flags[13] = malloc(4*sizeof(char));
     files[13] = get_string_item_(keywords, items, "RPSNAME_ROOT", nlines);
     if(files[13][0] == '\0'){ //no file
          strcpy(flags[13], "NO ");
     }
     else{
          strcpy(flags[13], "YES");
     }

     if(   (flags[10][0] == 'Y') || (flags[11][0] == 'Y') 
        || (flags[12][0] == 'Y') || (flags[13][0] == 'Y') ){
          files[9] = get_string_item_(keywords, items, "PSDIR", nlines);
     }
     else{
          files[9] = malloc(1*sizeof(char));
          files[9][0] = '\0';
     }

     /* END LEVINSON MODIFICATION FOR POSTAGE STAMP FILE OUTPUTS */

     /* BEGIN LEVINSON MODIFICATION FOR EXTRA FITTING LOOP */
     /* Ask if extra fitting loop tune3_.dofinalfit (int)*/
     char* dofinalfit_yn = get_string_item_(keywords, items, "DOFINALFIT", nlines);
     if ((dofinalfit_yn[0] == 'N' ) || (dofinalfit_yn[0] == 'n')){
          tune3_.dofinalfit = 0; //dont do final fit
     }
     else{
          tune3_.dofinalfit = 1; //if keyword is not specified or not NO, then yes 
     }

     /* END LEVINSON'S MODIFICATIONS */
   
     /* Done with flags and files extraction and update. 
        assign to correct tuneable varibales */
     
//printf("here26\n");
     for (I = 15; I < NFF; I++){
          files[I] = malloc(1*sizeof(char));
          files[I][0] = '\0';
          flags[I] = malloc(1*sizeof(char));
          flags[I][0] = '\0';
     }
     tune16_.files = files;
     tune16_.flags = flags;

//printf("here28\n");
     /* now write the parameters to file */
     /* add blank comments to any added keyword and value set */
     for(I = nlines_init; I < nlines; I++){
          comments[I] = malloc(2*sizeof(char));
          comments[I][0] = ' ';
          comments[I][1] = '\0';
     }
     /* write all keywords, items, and comments to file specified by value */
     paramwrite_(keywords, items, comments, &nlines);
//     printf("Parameters written to output file\n");

//printf("here end\n");
     /* free keywords, items, and comments arrays */
     for(I = 0; I < nlines; I++){
          free(keywords[I]);
          free(items[I]);
          free(comments[I]);
     }
     free(keywords);
     free(items);
     free(comments);

}

