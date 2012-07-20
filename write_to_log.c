#include <stdlib.h>
#include <stdio.h>
#include "io_dophotc.h"
#include "write_to_log.h"

/* routines to open and close the logfile */

/* global log file accessable by all other c routines */
FILE *logfile;

void openlog_(char* logname)
{
     char* prompt;
     logfile = fopen(logname, "w");
     if (logfile == NULL){
          fprintf(stderr, "logfile doesn't exist\n") ;
          free(logname);
          prompt = "Enter logfile name";
          logfile = opensc_(prompt, "w");
     }
}

void closelog_()
{
     fclose(logfile) ;
}

