#ifndef TUNEUP_H_
#define TUNEUP_H_

/* dophot subroutine converted to  c void function 04-01-2012 */
/* initializes all the tuneable parameters from the param files */

void tuneup_(int command_line_pm, char* pm_file_name);
/* command_line_pm is 1 if the parameter modification file was given at the 
   command line, 0 if not */

#endif
