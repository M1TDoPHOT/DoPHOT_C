#ifndef PARINTERP_STRUCT_H_
#define PARINTERP_STRUCT_H_

// having grown weary of the mallocs
//parinterp is called so often an by so many routines that
//it shouldn't share memory.
struct{
     short int* parsiarray; //[2]
     int*       parintarray; //[tune4_.npar-3]
     float*     dummy; //[npmax]
}pararray_;

#endif
