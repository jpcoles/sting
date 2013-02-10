#ifndef IC_H
#define IC_H

#include "sting.h"

int set_units(struct env *env);
int ic_random_circular(struct env *env);
int ic_random_elliptic(struct env *env);
int ic_2_particle_simple(struct env *env);
int ic_4_particle_simple(struct env *env);
int ic_coupled_pairs(struct env *env);
int ic_tube(struct env *env);
int ic_2tubes(struct env *env);
int ic_jet2d(struct env *env);


#endif
