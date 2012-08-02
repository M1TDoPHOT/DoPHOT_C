#ifndef OIMDEV_STRUCT_H
#define OIMDEV_STRUCT_H

struct
{
     float**   dxomp;
     float**   dyomp;
     int**       omp;
}oimdev_;
//DXOMP[2*IHSIDE][2*IHSIDE]
//DYOMP[2*IHSIDE][2*IHSIDE]
//OMP[2*IHSIDE+1][2*IHSIDE+1]

#endif
