#ifndef ADDOBJ_STRUCT_H
#define ADDOBJ_STRUCT_H

//Contains the four numbers describing the subtraction
// boundaries for each object in the field needed by
// addstar to add and subtract object from the image

struct{
     short int** addarea;
}addobj_;
//ADDAREA[NSMAX][4]

#endif
