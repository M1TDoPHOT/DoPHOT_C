#ifndef MODEL_STRUCT_H
#define MODEL_STRUCT_H

// contains the model information for each object
// whether pseudogaussian (1)
// or 7+ parameter or empirical (0) default

struct{
     int* which_model;
     int* tested;
}model_;
//which_model[NSMAX]
//tested[NSMAX]

#endif
