#ifndef PARAMFILE_H_
#define PARAMFILE_H_

/* dophot subroutine converted to c void functoin 03-31-2012 */

/* This subroutine asks the user for the parameter files.  If the modified file
   does not exist, the default file is transferred directly to the output file.
   A program (MKDEFAULT) is supplied in the DoPHOT directory to generate a 
   valid default parameter file should one not exist already. */

void paramfile_(int command_line_pm, char* pm_file_name, char** defkeyword, char** defitem, char** defcomment, int* nldef);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* dophot subroutine converted to c void functoin 03-31-2012 */
/* Write current parameters to unit PARAM_OUT file */
/* mostly a wrapper for writeitemc_ and fprintf that first finds the file name 
   in the item list */
/* command_line_pm is 1 if the parameter modification file was given at the
   command line, 0 if not */

void paramwrite_(char** defkeyword, char** defitem, char** defcomment, int* nldef);

#endif
