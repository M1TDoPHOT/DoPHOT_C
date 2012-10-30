#include <stdlib.h>
#include "logh.h"
#include "tuneable.h"
#include "free_parking_struct.h"
#include "search_struct.h"
#include "cast_arr.h"
#include "mini_mathlib.h"
#include "guess.h"
#include "oblims.h"
#include "toobright.h"

/* dophot logical function converted to c int function 02-03-2012 */

int toobright_(int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, float* STARPAR)
{

     /* dereference pointers */
     int NFAST = *NFAST_ptr;
     int NSLOW = *NSLOW_ptr;

     /* rename common block variables */
     short int* KRECT = tune2_.krect ;
     int ITOP       = tune3_.itop    ; 
     int ICRIT      = tune6_.icrit   ;
     float CTPERSAT = tune6_.ctpersat;
     float CMAX     = tune6_.cmax    ;
     int N0LEFT     = tune7_.n0left  ;
     int N0RIGHT    = tune7_.n0right ;
     int lverb      = tune14_.lverb  ;


     /* substance of function begins here */
     float* A     = free_parking_.npmaxarray_1;
     short* JRECT = free_parking_.sifourarray_1;
     int FLAG;
     int JTOP, NSAT; 
     int IXHI, IXLO, IYHI, IYLO;
     int IX, IY, JY, JX;
     double DUM;
     int TOOBRIGHT = 0; //false

     if (STARPAR[1] < (float)(ITOP)/4.0f){
          return TOOBRIGHT;
     }
     else{
//          JTOP = min(ITOP, 32767);
          JTOP = min(ITOP, MAGIC);
          NSAT = 0;
          DUM  = guess2_(A, STARPAR, &IX, &IY);
          IXHI = min((int)(IX + KRECT[0]/2 + 0.5f), NFAST - N0RIGHT);
          IXLO = max((int)(IX - KRECT[0]/2 + 0.5f), 1     + N0LEFT );
          IYHI = min((int)(IY + KRECT[1]/2 + 0.5f), NSLOW);
          IYLO = max((int)(IY - KRECT[1]/2 + 0.5f), 1    );
          for (JY = (IYLO - 1); JY < IYHI; JY++){
               for (JX = (IXLO - 1); JX < IXHI; JX++){
                    if ((NOISE[JY][JX] == MAGIC) &&
                        (  BIG[JY][JX] >= JTOP )) {
                         NSAT += 1;
                    }
               }
          }
     }
     /* flag if number of saturated pixels if NSAT is too large */
     FLAG = ((STARPAR[1] > CMAX) || (NSAT >= ICRIT));
     if (FLAG){
          TOOBRIGHT = 1; //TRUE
          if (NSAT >= ICRIT){
               STARPAR[1] = CTPERSAT*NSAT;
               if (lverb > 10){
                    fprintf(logfile,"SATURATED PIXELS in object at %d %d\n",
                            IX, IY);
               }
          }
          else{
               STARPAR[1] = CTPERSAT*0.25;
          }
          oblims_(STARPAR, JRECT);
          STARPAR[4] = JRECT[1] - JRECT[0];
          STARPAR[5] = 0.01f;  //0 exactly saved for cosmic rays
          STARPAR[6] = JRECT[3] - JRECT[2];
     }

     return TOOBRIGHT;
}
 

