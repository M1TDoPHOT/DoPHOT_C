#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tuneable.h"
#include "io_dophotc.h"
#include "cast_arr.h"
#include "paramfile.h"

/* dophot subroutine converted to c void functoin 03-31-2012 */

/* This subroutine asks the user for the parameter files.  
   and returns the full set of default parameters modified where
   necessary by the modified parameters file. The modified params
   file is requested and the default params file is found in
   the mod param file or requested if not found */
/* keyword, item, and comment lists contents are malloced here,
   byt the list of pointers must be malloced elsewhere */

void paramfile_(int command_line_pm, char* pm_file_name, char** defkeyword, char** defitem, char** defcomment, int* nldef)
{
     //nphead defined to 400 in tuneup.c for consistency

     char headerline[180];
     char* keyword; //malloc'd to proper size in readitemc 
     char* item   ; //malloc'd to proper size in readitemc
     char* comment; //malloc'd to proper size in readitemc
     //from mod_param_file
     FILE*  mod_paramfile; //modified parameter file
     //individual strs are malloced seperately to correct sizes
     char** modkeyword = malloc(nphead*sizeof(char*)); 
     char** moditem    = malloc(nphead*sizeof(char*));
     char** modcomment = malloc(nphead*sizeof(char*));
     int nlmod;
     //from default_param_file
     FILE*  def_paramfile; //default parameter file
     char*  def_paramfile_name; //malloc'd to proper size in stringstrip

     char* prompt;
     char* rwa  = "r";
     int notfound, readerr, goodline;
     int still_needs_file;
     int I, J;

     /* ask for an open the mod_paramsfile by calling opens */
     if (command_line_pm == 1){ /* parameter modifications file given 
                                   at command linea */
          mod_paramfile = openac_(&readerr, pm_file_name, rwa);
          if (readerr != 0){
               printf("Parameter Modification file input on command line ");
               printf("not found!\n");
          }
     }
     if ((command_line_pm == 0) || (readerr == 1)){
          readerr = 0;
          prompt = "Enter modification parameters file name";
          mod_paramfile = opensc_(prompt, rwa);
     }
          

     /* collect lines from mod_param_file */ 
     I = 0;
     while( fgets(headerline, 180, mod_paramfile) != NULL ){
          //readitem mallocs keyword, item, comment
          readitemc_(headerline, &keyword, &item, &comment, &goodline);
          if (goodline){ //keyword exists
               modkeyword[I] = keyword;
               moditem[I]    = item;
               modcomment[I] = comment;
               I++;
          }
          else{
               free(keyword);
               free(item);
               free(comment);
          }
     }
     nlmod = I;
     
     /* find out if the names of the default and output parameter
        files can be found in mod_param_file.  If not, ask for these. */
     I = 0;
     notfound = 1;
     while ((notfound) && (I < nlmod)){
          if (strcmp(modkeyword[I], "PARAMS_DEFAULT") == 0){ 
               notfound = 0;
               readerr = 0; //default assume no error
               if (strlen(moditem[I]) > 0){ //item assigned to keyword
                    //strip name of leading spaces and quotes
                    def_paramfile_name = stringstrip_(moditem[I]);
                    //open file of given name
                    def_paramfile = openac_(&readerr, def_paramfile_name, rwa);
                    free(def_paramfile_name);
               }
               if ((strlen(moditem[I]) <= 0) || (readerr == 1)){
                    free(moditem[I]);
                    printf("Default parameters file not found!\n");
                    prompt = "Enter default parameters file name";
                    //NOTE:  I'm essentially rewritting opensc here because
                    //       I want to keep the image name
                    still_needs_file = 1; //true
                    while (still_needs_file){
                         printf("%s : ", prompt);
                         if (scanf("%s", def_paramfile_name) != 1){
                              printf("Error in file name.  Try again.\n");
                         }
                         else{
                              def_paramfile = fopen(def_paramfile_name, rwa);
                              if (def_paramfile == NULL){
                                   printf("Error opening file.  Try again.\n");
                              }
                              else{
                                   still_needs_file = 0;
                              }
                         }
                    }
                    moditem[I] =
                         malloc((strlen(def_paramfile_name)+3)*sizeof(char));
                    sprintf(moditem[I], "\'%s\'", def_paramfile_name); 
               }
          }
          I++;
     }
     if (notfound){
          def_paramfile_name = malloc(64*sizeof(char));
          printf("Default parameters file not found!\n");
          prompt = "Enter default parameters file name";
          def_paramfile = opensc_(prompt, rwa);
     }

     /* Loop over the default parameter file and look for those keywords 
        in the modified parameter file. */
     I = 0;
     /* collect lines from mod_param_file */ 
     while( fgets(headerline, 180, def_paramfile) != NULL){
          readitemc_(headerline, &keyword, &item, &comment, &goodline);
          if (goodline){
               defkeyword[I] = keyword;
               /* See if you can find the keyword and val in the modified 
                  parameter file. If not write the default header variables
                  without change */
               notfound = 1; //true
               J = 0;
               while ((notfound) && (J < nlmod-1)){
                 if (strcmp(keyword, modkeyword[J]) == 0){
                   notfound = 0;
                   if (strlen(moditem[J]) > 0){
                     defitem[I] = malloc((strlen(moditem[J])+1)*sizeof(char));
                     strcpy(defitem[I], moditem[J]);
                     free(item);
                   }
                   else{
                     defitem[I] = item;
                   }
                   if (strlen(modcomment[J]) > 1){
                     defcomment[I] = malloc((strlen(modcomment[J])+1)*sizeof(char));
                     strcpy(defcomment[I], modcomment[J]);
                     free(comment);
                   }
                   else{
                     defcomment[I] = comment;
                   }
                 }
                 J++;
               } //end loop over mod file keywords
               if (notfound){ // if keyword not found, write default vals
                    defitem[I] = item;
                    defcomment[I] = comment;
               }
               I++;
          }
          else{ //if bad line, free allocated memory
               free(keyword);
               free(item);
               free(comment);
          }
     }
     *nldef = I;
     /* free locally allocated memory */
     fclose(mod_paramfile);
     fclose(def_paramfile);

     free_char_arr(nlmod, modkeyword);
     free_char_arr(nlmod, moditem);
     free_char_arr(nlmod, modcomment);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* dophot subroutine converted to c void functoin 03-31-2012 */
/* Write current parameters to unit PARAM_OUT file */
/* mostly a wrapper for writeitemc_ and fprintf that first finds the file name 
   in the item list */

void paramwrite_(char** defkeyword, char** defitem, char** defcomment, int* nldef)
{

     FILE*  out_paramfile; //modified parameter file
     char*  out_paramfile_name;
     char*  prompt;
     char line[180];
     int I;
     int notfound, readerr, still_needs_file;
     char* rwa = "w";

     /* get output parameters file from defitems and open it
        If error, ask for new file */
     I = 0;
     notfound = 1;
     while( (notfound) && (I < *nldef) ){
          if(strcmp(defkeyword[I], "PARAMS_OUT") == 0){
               notfound = 0;
               readerr = 0; //default assume no error
               if (strlen(defitem[I]) > 0){ //item assigned to keyword
                    //strip name of leading spaces and quotes
                    out_paramfile_name = stringstrip_(defitem[I]);
                    //open file of given name
                    out_paramfile = openac_(&readerr, out_paramfile_name, rwa);
                    free(out_paramfile_name);
               }
               if ((strlen(defitem[I]) <= 0) || (readerr == 1)){
                    free(defitem[I]);
                    printf("Error opening output parameters file.\n");
                    prompt = "Enter output parameters file name";
                    //NOTE:  I'm essentially rewritting opensc here because
                    //       I want to keep the image name
                    still_needs_file = 1; //true
                    while (still_needs_file){
                         printf("%s : ", prompt);
                         if (scanf("%s", out_paramfile_name) != 1){
                              printf("Error in file name.  Try again.\n");
                         }
                         else{
                              out_paramfile = fopen(out_paramfile_name, rwa);
                              if (out_paramfile == NULL){
                                   printf("Error opening file.  Try again.\n");
                              }
                              else{
                                   still_needs_file = 0;
                              }
                         }
                    }
                    defitem[I] =
                         malloc((strlen(out_paramfile_name)+3)*sizeof(char));
                    sprintf(defitem[I], "\'%s\'", out_paramfile_name); 
               }
          }
          I++;
     }

     // write keyword, val, comment to a single line and then to lines in file
     for(I = 0; I < *nldef; I++){
          writeitemc_(line, defkeyword[I], defitem[I], defcomment[I]);
          fputs(line, out_paramfile);
     }

     //close file, free locally allocated memory
     fclose(out_paramfile);
}

 
