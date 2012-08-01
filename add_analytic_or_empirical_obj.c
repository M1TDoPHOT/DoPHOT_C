#include <stdlib.h>
#include "tuneable.h"
#include "estuff_struct.h"
#include "galpass_struct.h"
#include "addstar.h"
#include "empiricals.h" //contains oldemp and oneemp
#include "add_analytic_or_empirical_obj.h"

/* addstar wrapper: switch between model and footprint to add using addstar */

void add_analytic_or_empirical_obj(double (*analytic_model)(short int*, float*, float*, int*, int*), int** big, int** noise, int nfast, int nslow, float** starpar, short int** addarea, int add_or_subtract, int model_img, char* model_file, int clean_img, char* clean_file, int obj_index, int galaxy)
{

     /* add_or_subtract is 1 (any positive) for add and 
        -1 (any negative) for subtract */
     /* galaxy is 1 for a galaxy and 0 else */

     int iadd; //passed to addstar
     float* pass_par = starpar[obj_index];
     short int* pass_area = addarea[obj_index];
     double (*pass_model)(short int*, float*, float*, int*, int*) = analytic_model;

     if ( (tune22_.emenab) && (tune22_.empok) && 
          (estuff_.emsub[obj_index] == 1) ){
          pass_par   = estuff_.empar[obj_index];
          pass_model = &oneemp_;

          if (add_or_subtract > 0) iadd = 2;
          else if (add_or_subtract < 0) iadd = -2;

     }
     else{ //analytic or empirical template

          if (add_or_subtract > 0) iadd = 1;
          else if (add_or_subtract < 0) iadd = -1;

          if (galaxy){
               galpass_.bigfoot = 1;
          }

     }

     addstar_(pass_model, big, noise,
                   nfast, nslow,
                   pass_par, pass_area, iadd,
                   model_img, model_file, 
                   clean_img, clean_file); 

     galpass_.bigfoot = 0; //false

}
