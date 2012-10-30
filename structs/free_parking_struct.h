#ifndef FREE_PARKING_STRUCT_H_
#define FREE_PARKING_STRUCT_H_

// having grown weary of the mallocs, here is a list of
// arrays which are to only live inside of individual routines
// and NEVER have values saved to them.  all needed values 
// should be set and not assumed to be 0.
struct{
     float* npmaxarray_1;
     float* npmaxarray_2;
     float* npmaxarray_3;
     float* npmaxarray_4;
     short int* sifourarray_1;
     short int* sifourarray_2;
     int* intfourarray_1;
     int* intfourarray_2;
     float** npmaxbynpmaxarray;
}free_parking_;

#endif
