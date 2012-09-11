#ifndef STARLIST_STRUCT_H
#define STARLIST_STRUCT_H

struct{
     float** starpar;
     int* IMTYPE;
     float** shadow;
     float** shaderr;
     float*** shadcovar;
}starlist_;
//STARPAR[NSMAX][NPMAX]
//SHADOW[NSMAX][NPMAX]
//SHADERR[NSMAX][NPMAX]
//SHADCOVAR[NSMAX][NPMAX][NPMAX]

#endif
