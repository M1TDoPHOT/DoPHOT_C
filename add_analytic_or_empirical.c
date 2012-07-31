#include <stdlib.h>
#include "tuneable.h"
#include "addstar.h"
#include "empiricals.h" //contains oldemp and oneemp
#include "add_analytic_or_empirical_obj.h"

/* addstar wrapper: switch between model and footprint to add using addstar */

void add_analytic_or_empirical_obj(double (*MODEL)(short int*, float*, float*, int*, int*), int** BIG, int** NOISE, int* NFAST_ptr, int* NSLOW_ptr, float** STARPAR, short int** ADDAREA, int ADD, int model_img, char* model_file, int clean_img, char* clean_file, int EMSUB)
{

     /* no need to dereferenec pointer or rename used common block vars */
     if ((tune22_.emenab) & (tune22_.empok)){
          if (EMSUB >= 1){
               addstar_(&oneemp_, BIG, NOISE,
                        &NFAST, &NSLOW,
                        EMPAR[K],
                        ADDAREA[K], &JADD,
                        0, " ", 0, " ");
     }
     else{
          galpass_.bigfoot = (JMTYPE == 2);
          addstar_(ONESTAR, BIG, NOISE,
                   &NFAST, &NSLOW,
                   STARPAR[K],
                   ADDAREA[K], &IADD,
                   0, " ", 0, " ");
          galpass_.bigfoot = 0; //false
     }


