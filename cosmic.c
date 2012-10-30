#include <stdlib.h>
#include "logh.h"
#include "tuneable.h"
#include "free_parking_struct.h"
#include "fitarrays_struct.h"
#include "unitize_struct.h"
#include "mini_mathlib.h"
#include "cast_arr.h"
#include "guess.h"
#include "cosmic.h"

/* dophot logical function converted to c int function 02-21-2012 */
//onestar is a pointer to the function name ONESTAR which returns double
//onestar is usually pseud2d
int cosmic_( double (*ONESTAR)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, float* STARPAR )
{

     /* dereference pointers */
     int NFAST = *NFAST_ptr;
     int NSLOW = *NSLOW_ptr;
     /* rename used common block variables */
     int NPAR = tune4_.npar;
     float WIDOBL  = tune6_.widobl ;
     float DISCRIM = tune8_.discrim;
     float SN2COS  = tune9_.sn2cos ;
     int lverb = tune14_.lverb;
     float* A  = fitarrays_.a ;
     float* FA = fitarrays_.fa;
     float UFACTOR = unitize_.ufactor;

     /* substance of function begins here */
     short int* II = free_parking_.sifourarray_1;
     double guess2_return;
     double onestar_return;
     float SKY, CHISTAR, CHICOS;
     float MAXOB;
     int NPIX;
     int IX, IY;
     int KX, KY, condition;
     float PRED, OBS, TEMP, SN2;
     int IMAX, JMAX;
     float TNOISE;
     int POINTY = 0; //false
     int COSMIC = 0; //false 
     int ZERO_dum = 0; //allows 0 to be passed as pointer

     guess2_return = guess2_(A, STARPAR, &IX, &IY);
     SKY     = (float)guess2_return;
// short int data remnant     MAXOB   = -32768.0f;
     MAXOB   = (float)(OBLITVAL);
     CHISTAR = 0.0f;
     CHICOS  = 0.0f;
     NPIX    = 0;
     for(KY = -1; KY <= 1; KY++){
          for(KX = -1; KX <= 1; KX++){
               condition = (  ((IX+KX) <= NFAST)
                           && ((IX+KX) >= 1    )
                           && ((IY+KY) <= NSLOW)
                           && ((IY+KY) >= 1    ) );
               if (condition){
                    II[0] = KX;
                    II[1] = KY;
                    if (NOISE[IY+KY-1][IX+KX-1] < MAGIC){
                         NPIX += 1;
                         onestar_return =(*ONESTAR)(II,A,FA,&NPAR,&ZERO_dum);
                         PRED = UFACTOR*((float)(onestar_return));
                         OBS  = (float)(BIG[IY+KY-1][IX+KX-1]);
                         TEMP = 1.0f/((float)(NOISE[IY+KY-1][IX+KX-1]));
                         CHISTAR += (OBS - PRED)*(OBS - PRED)*TEMP;
                         SN2 = (OBS - SKY)*(OBS - SKY)*TEMP;
                         CHICOS += SN2;
                         if ( (OBS > MAXOB) && (SN2 >= SN2COS) ){
                              IMAX  = IX + KX; 
                              JMAX  = IY + KY; 
                              MAXOB = OBS; 
                              TNOISE = TEMP;
                         }
                    }
               } // end of 'condition' if
          }//end kx loop
     }// end ky loop
                         
//     if ((NPIX > 7) && (MAXOB != -32768.0f)){
     if ((NPIX > 7) && (MAXOB != (float)(OBLITVAL))){
          CHICOS -= (MAXOB - SKY)*(MAXOB - SKY)*TNOISE;
          POINTY  = (CHICOS < DISCRIM*CHISTAR);
          if (lverb > 30){
               fprintf(logfile,"Location: %d %d\n", IX, IY);
               fprintf(logfile,"CHI-STAR & CHI-COSMIC = %f %f\n", CHISTAR, CHICOS);
          }
     }
     if (POINTY){
          if (lverb > 20){
               fprintf(logfile,"COSMIC RAY INTENSITY, IX, IY : %d %d %d\n",
                    BIG[JMAX-1][IMAX-1], IX, IY);
          }
          STARPAR[0] = SKY         ;
          STARPAR[1] = MAXOB       ;
          STARPAR[2] = (float)IMAX ;
          STARPAR[3] = (float)JMAX ;
          STARPAR[4] = WIDOBL      ;
          STARPAR[5] = -1.0f       ;
          STARPAR[6] = WIDOBL      ;
     }
     COSMIC = POINTY;
     
     return COSMIC;
}
     

