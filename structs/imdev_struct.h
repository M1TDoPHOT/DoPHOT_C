#ifndef IMDEV_STRUCT_H
#define IMDEV_STRUCT_H

struct
{
     float**   dxemp;
     float**   dyemp;
     int**       emp;
}imdev_;
//DXEMP[2*IHSIDE][2*IHSIDE]
//DYEMP[2*IHSIDE][2*IHSIDE]
//EMP[2*IHSIDE+1][2*IHSIDE+1]

#endif
