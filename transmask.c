#include "logh.h"
#include "tuneable.h"
#include "passmask_struct.h"
#include "trans7_struct.h"
#include "passimp_struct.h"
#include "cast_arr.h"
#include "transmask.h"

int transmask_(int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, int* IX_ptr, int* IY_ptr, float* SKY_ptr, float* HUMP_ptr)
{

     /* dereference pointers */
     int NFAST  = *NFAST_ptr;
     int NSLOW  = *NSLOW_ptr;
     int IX     = *IX_ptr   ;
     int IY     = *IY_ptr   ;
     float SKY  = *SKY_ptr  ;
//unused     float HUMP = *HUMP_ptr ;

     /* rename used common block variables */
     int   TEST7 = trans7_.test7 ;
     int   III   = passimp_.iii  ;
     int   IXBY2 = tune12_.ixby2 ;
     int   IYBY2 = tune12_.iyby2 ;
     int   lverb = tune14_.lverb ;
     float  CRIT7     = tune9_.crit7    ;
     float  BUMPCRIT  = tune9_.bumpcrit ;
     float** STARMASK = passmask_.starmask;
     
     /* substance of function begins here */
     float HUMP2, PMET;
     float SUM0, SUM1, SUM2;
     float RAT;
     int I, J, II, JJ;
     int TRANSMASK = 0; //false

     int this_big;
     int this_noise;
     float this_starmask;

     if (TEST7){
          HUMP2 = CRIT7;
     }
     else{
          HUMP2 = BUMPCRIT*BUMPCRIT;
     }
     SUM0 = 0.0f;
     SUM1 = 0.0f;
     SUM2 = 0.0f;
     // loop is still fortran style 1 indexed for ease of conversion
     // so all array references are i-1 or j-1
     for (J = -IYBY2; J <= IYBY2; J++){
          JJ = IY + J;
          if ( (JJ >= 1) && (JJ <= NSLOW) ){
               for (I = -IXBY2; I <= IXBY2; I++){
                    II = IX + I;
                    if ( (II >= 1) && (II <= NFAST) ){ 
                         this_noise    = NOISE[JJ-1][II-1];
                         this_starmask = STARMASK[J+IYBY2][I+IXBY2]; 
                         this_big      = BIG[JJ-1][II-1];
                         if (this_noise != MAGIC){
                              PMET = this_starmask/(float)this_noise;
                              SUM0 += 1.0f;
                              SUM1 += PMET*this_starmask;
                              SUM2 += PMET*( (float)(this_big) - SKY );
                         }
                    }
               }// end I loop
          }
     }// end J loop

     if (SUM2 > 0.0f){
          RAT = SUM2*SUM2 / SUM1;
          if (RAT > HUMP2){
               TRANSMASK = 1; //true
               if ( RAT < (1.1f*HUMP2) ){
                    if (lverb > 30) {
                         fprintf(logfile, "MARGINAL: (S/N)**2 through Mask = %f \n",
                                      RAT);
                    }
               }
          }
     }
     else{
          if (III == 7){
               fprintf(logfile, "RISKY: NO GOOD PIXELS FOR TRANSMASK:\n");
               fprintf(logfile, "IX, IY = %d %d \n", IX, IY);
               fprintf(logfile, "sum0, sum1, sum2 = %f %f %f \n", 
                                 SUM0, SUM1, SUM2);
          }
     }

     return TRANSMASK;
}
