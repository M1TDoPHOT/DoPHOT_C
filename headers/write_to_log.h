#ifndef WRITE_TO_LOG_H_
#define WRITE_TO_LOG_H_

/* routines to open and close the logfile */

/* global log file accessable by all other c routines */
FILE *logfile;

void openlog_(char* logname);

void closelog_();

#endif
