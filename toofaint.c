#include "logh.h"
#include "tuneable.h"
#include "unitize_struct.h"
#include "toofaint.h"

/* dophot function converted to c int function 02-02-2012 */
/* determines if star is too faint for analysis */

int toofaint_(float* STARPAR, float* ERR)
{
    
     /* no need to dereference pointers as already 1d arrays */
     /* simplifying common block variable names */
     float UFACTOR = unitize_.ufactor;
     float CRIT7 = tune9_.crit7;
     int lverb = tune14_.lverb;

     /* substance of function begins here */
     int ret_val ;
     float SIG2 ;

     ret_val = (STARPAR[1] < 1.0f);
     if (ERR[1] > 0){
          SIG2 = (STARPAR[1]/UFACTOR)*(STARPAR[1]/UFACTOR) / ERR[1] ;
          ret_val = ( ret_val || (SIG2 < CRIT7) ) ;
          if (lverb > 30){
               fprintf(logfile,"SIG2 = %9.6f\n", SIG2);
          }
     }

     return ret_val;
}

	
