#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cast_arr.h"
#include "io_dophotc.h"

/* subroutines to open and close files, get header info, etc */
/* these routines will not interact with fortran, but rather are only
   callable by c programs */

/*
c  This subroutine opens or closes a sequential file
c    for reading by taking a file name at the terminal.
c  It also handles basic errors.
// replaces fortran 'opens'
*/
FILE* opensc_(char* prompt, char* rwa)
{
     FILE* this_file;
     char file_name[80];
     int still_needs_file = 1; //true
     while (still_needs_file){
          printf("%s : ", prompt);
          if (scanf("%s", file_name) != 1){
               printf("Error in file name.  Try again.\n");
          }
          else{
               this_file = fopen(file_name, rwa);
               if (this_file == NULL){
                    printf("Error opening file.  Try again.\n");
               }
               else{
                    still_needs_file = 0; //false, exit routine
               }
          }
     }
     return this_file;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*
c  This subroutine opens or closes a sequential file given a file name.
c  It returns an error code, 0 if ok, 1 if there is a problem.
// replaces fortran 'opena'
*/
FILE* openac_(int* ierr, char* file_name, char* rwa)
{
     *ierr = 0; //ok
     FILE* this_file = fopen(file_name, rwa);
     if (this_file == NULL){
          *ierr = 1; //not ok
     }
     return this_file;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*
c  Header reading/writing subroutines.
*/

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* clean string */
void charinitc_(char* str, int len)
{
     memset(str, '\0', len);
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*
c  Reads a line and extracts the keyword, parameter, and comment if they exist
c  from a single header line.  also returns the item type as an int.
*/
// combines the functionality of the fortran version's findkey, getcomment,
// and readitem
// WARNING mallocs!!!!
// much of the logic credit here goes to JAL who was bored and considered this
// "real programming"
void readitemc_(char* line, char** keyword, char** value, char** comment, int* goodline)
{

  char *begin;
  char *end;

  //only set any values not null if
  //and doesn't start with an equals (no keyword, comment line)
  if ((strlen(line) >= 2) && (line[0] != '=')){ 
    *goodline = 1;  
    end = line;
  
    //mark beginning of end of leading spaces
    while (*end == ' ') ++end;
    begin = end;
    //mark end of first word
    while (*end != ' ' && *end != '=' && *end != '\0') ++end;

    //set first word to *keyword
    *keyword = malloc(sizeof(char) * (end - begin + 1));
    strncpy(*keyword, begin, end - begin);
    keyword[0][end - begin] = '\0';

    //mark beginning of end of next set of leading spaces
    while (*end == ' ') ++end;

    //if equals sign found, find end of spaced after and set to value
    //else set *value equal to null character
    if (*end == '=') {
      do {
        ++end;
      } while (*end == ' ');
      begin = end;
      while (*end != ' ' && *end != '\0') ++end;

      *value = malloc(sizeof(char) * (end - begin + 1));
      strncpy(*value, begin, end - begin);
      value[0][end - begin] = '\0';
    }
    else {
      *value  = malloc(sizeof(char) * 1);
      **value = '\0';
    }

    //mark beginning of end of next set of leading spaces
    while (*end == ' ') ++end;
    //grab everything else until ending spaces as comment
    //if it exists
    if (*end != '\0') {
    begin = end;
      while (*end) ++end;
      while (*(end - 1) == ' ') --end;

      *comment = malloc(sizeof(char) * (end - begin + 1));
      strncpy(*comment, begin, end - begin);
      comment[0][end - begin] = '\0';
    }
    else {
      *comment = malloc(sizeof(char) * 1);
      **comment = '\0';
    }

  } 
  else{ //if comment line or blank line
    *goodline = 0; 
    *keyword = malloc(sizeof(char) * 1);
    *value   = malloc(sizeof(char) * 1);
    *comment = malloc(sizeof(char) * 1);
    **keyword = '\0';
    **value   = '\0';
    **comment = '\0';
  }

//  printf("keyword : %s\n", *keyword);
//  printf("moditem : %s\n", *value);
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*
   determines an item type from its char nature and its keyword (int v. float)
*/
/* itype = 0: float  */
/* itype = 1: int    */
/* itype = 2: string */
void typeitem_(char* keyword, char* value, int* itype)
{
     char first;
     //determine if string variable by examining value for leading '
     if (value[0] == '\''){
          *itype = 2;
     }
     else{ // examine keyword name for indicators of int v. float
          first = keyword[0];
          if ((first == 'I') || (first == 'i') ||
              (first == 'J') || (first == 'j') ||
              (first == 'K') || (first == 'k') ||
              (first == 'L') || (first == 'l') ||
              (first == 'M') || (first == 'm') ||
              (first == 'N') || (first == 'n') ){
               *itype = 1;
          }
          else{ 
               *itype = 0;
          }
     }  
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* strip a string item of leading spaces and quotes */
//MALLOCS
char* stringstrip_(char* instr)
{
     char* begin;
     char* end;
     end = instr;
     //strip of leading spaces and quotes
     while((*end == '\'') || (*end == ' ')) ++end;
     begin = end;
     //mark the end given by quotes or space or null
     while((*end != '\'') && (*end != ' ') && (*end != '\0')) ++end;
     //set outstring to that inbetween markings
     char* outstr = malloc(sizeof(char) * (end - begin + 1));
     strncpy(outstr, begin, end - begin);
     outstr[end - begin] = '\0';
     return outstr;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* get a string item from the list of keywords, stripping it of quotes */
// MALLOCS via stringstrip
char* get_string_item_(char** keywords, char** items, char* keyword, int len)
{
     int I = 0;
     int notfound = 1;
     char* item;
     while ((notfound) && (I < len)){
          if (strcmp(keywords[I], keyword) == 0){
               notfound = 0;
               if (items[I][0] != '\0'){
                    item = stringstrip_(items[I]);
               }
               else{
                    item = malloc(1*sizeof(char));
                    *item = '\0';
               }
          }
          I++;
     }
     if (notfound){
          item = malloc(1*sizeof(char));
          *item = '\0';
     }
     return item;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* put a string item into the correct place given the list of keywords,
   if type = 2, assume string and add quotes */
// if keyword doesn't exits, add 1 to len and add keyword and item
void put_item_(char** keywords, char** items, char* keyword, int *len, char* newitem, int type)
{
     int I = 0;
     int notfound = 1;
     char* temp = malloc(64*sizeof(char));
     while ((notfound) && (I < *len)){
          if (strcmp(keywords[I], keyword) == 0){
               notfound = 0;
               free(items[I]);
               if(type == 2){ //for string type, add quotes
                    sscanf(temp, "\'%s\'",newitem);
                    items[I] = malloc((strlen(newitem)+3)*sizeof(char));
                    strcpy(items[I], temp);
               }
               else{
                    items[I] = malloc((strlen(newitem)+1)*sizeof(char));
                    strcpy(items[I], newitem);
               }
          }
          I++;
     }
     if (notfound){
          keywords[I] = malloc(strlen(keyword)*sizeof(char));
          strcpy(keywords[I], keyword);
          if(type == 2){ //for string type, add quotes
               sscanf(temp, "\'%s\'", newitem);
               items[I] = malloc((strlen(newitem)+3)*sizeof(char));
               strcpy(items[I], temp);
          }
          else{
               items[I] = malloc((strlen(newitem)+1)*sizeof(char));
               strcpy(items[I], newitem);
          }
          *len += 1;         
     }
     free(temp);
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* get an int item from the list of keywords */
int get_int_item_(char** keywords, char** items, char* keyword, int len)
{
     int I = 0;
     int notfound = 1;
     int item;
     while ((notfound) && (I < len)){
          if (strcmp(keywords[I], keyword) == 0){
               notfound = 0;
               if (items[I][0] != '\0'){
                    item = atoi(items[I]); 
               }
               else{
                    item = -999;
               }
          }
          I++;
     }
     if (notfound){
          item = -999;
     }
     return item;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* get a float item from the list of keywords */
float get_float_item_(char** keywords, char** items, char* keyword, int len)
{
     int I = 0;
     int notfound = 1;
     float item;
     while ((notfound) && (I < len)){
          if (strcmp(keywords[I], keyword) == 0){
               notfound = 0;
               if (items[I][0] != '\0'){
                    item = atof(items[I]); 
               }
               else{
                    item = -999;
               }
          }
          I++;
     }
     if (notfound){
          item = -999.9f;
     }
     return item;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*
c  Writes values out to a header line.  Also writes the comment for that
c  line if ncomm .gt. 0.  nr is the length of the item string.  ITYPE is
c  coded as described in RFITEM.
c  REQUIRES knowing passing item as string type, so need to convert
c  item to string before call, can be done with sprintf
*/
void writeitemc_(char* line, char* keyword, char* value, char* comment)
{
     
  strcpy(line, keyword);
  if (strlen(keyword) > 0 ){
    strcat(line, " = ");
    strcat(line, value);
  }
  if (strlen(comment) > 1 ){
    strcat(line, "          "); //inserting 10 spaced before comment
    strcat(line, comment);
  }

}
  
