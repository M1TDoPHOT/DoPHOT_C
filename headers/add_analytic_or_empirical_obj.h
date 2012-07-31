#ifndef ADD_ANALYTIC_OR_EMPIRICAL
#define ADD_ANALYTIC_OR_EMPIRICAL

/* addstar wrapper: switch between model and footprint to add using addstar */

void add_analytic_or_empirical_obj(double (*analytic_model)(short int*, float*, float*, int*, int*), int** big, int** noise, int nfast, int nslow, float** starpar, short int** addarea, int add_or_subtract, int model_img, char* model_file, int clean_img, char* clean_file, int obj_index, int galaxy);

#endif
