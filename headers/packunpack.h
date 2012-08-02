#ifndef PACKUNPACK_H_
#define PACKUNPACK_H_

/* this abomination is designed to take two shorts and pack them into
   one float and then pull the two shorts back out.  it is a left over
   from fortran's equivalence statement.  It cannot be killed until 
   all routines that use it 
	addstar, chisq, cosmic, fillerup, probgal, and varipar
   no longer do so.  moreover, the floats it converts shorts into are
   often fed into an array, XX in the struct 'subraster' which must also
   be converted to be a 2d array of ints rather than a 1d array of floats.
   and god help us all if anything is equated to subraster_.xx...
   boys and girls, this is why we dont use globals */

void pack_(short int* jx, short int* jy, float* r4);
    
void unpack_(short int* jx, short int* jy, float* r4);

long mask_();
   
#endif 
