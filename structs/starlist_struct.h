#ifndef STARLIST_STRUCT_H
#define STARLIST_STRUCT_H

struct{
     float** starpar;
     int* IMTYPE;
     float** shadow;
     float** shaderr;
}starlist_;
//STARPAR[NSMAX][NPMAX]
//SHADOW[NSMAX][NPMAX]
//SHADERR[NSMAX][NPMAX]

#endif
