#ifndef IO_DOPHOTC_H_
#define IO_DOPHOTC_H_

/* subroutines to open and close files, get header info, etc */
/* loosely translated from fortran, but mostly rewritten for c */
/* these routines will not interact with fortran, but rather are only
   callable by c programs */

/*
c  This subroutine opens or closes a sequential file
c    for reading by taking a file name at the terminal.
c  It also handles basic errors.
// replaces fortran 'opens'
*/
FILE* opensc_(char* prompt, char* rwa);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*
c  This subroutine opens or closes a sequential file given a file name.
c  It returns an error code, 0 if ok, 1 if there is a problem.
// replaces fortran 'opena'
*/
FILE* openac_(int* ierr, char* file_name, char* rwa);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* clean string */
void charinitc_(char* str, int len);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*
c  Reads a line and extracts the keyword, item/value, and comment if it exists
c  from a single header line.  
*/
//mallocs keyword, value, and comment within the routine to the proper length
//c handles strings gracelessly.  keyword is a pointer to a string (char*)

void readitemc_(char* header, char** keyword, char** value, char** comment, int* goodline);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*
   determines an item type from its char nature and its keyword (int v. float)
*/
/* itype = 0: float  */
/* itype = 1: int    */
/* itype = 2: string */
void typeitem_(char* keyword, char* value, int* itype);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* strip a string item of leading spaces and quotes */
//mallocs its return
char* stringstrip_(char* instr);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* get a string item from the list of keywords, stripping it of quotes */
// MALLOCS via stringstrip
char* get_string_item_(char** keywords, char** items, char* keyword, int len);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* get an int  item from the list of keywords */
int get_int_item_(char** keywords, char** items, char* keyword, int len);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* get a float item from the list of keywords */
float get_float_item_(char** keywords, char** items, char* keyword, int len);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* put a item into the correct place given the list of keywords,
   adding the quotes for a string item, type = 2*/
// if keyword doesn't exits, add 1 to len and add keyword and item
void put_item_(char** keywords, char** items, char* keyword, int* len, char* newitem, int type);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*
c  Writes values out to a header line.  Also writes the comment for that
c  line if ncomm .gt. 0.  
c  REQUIRES knowing passing item as string type, so need to convert
c  item to string before call 
*/
void writeitemc_(char* line, char* keyword, char* value, char* comment);

#endif
