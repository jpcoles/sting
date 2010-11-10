#ifndef STING_H
#define STING_H

#include <stdio.h>

#define WITH_FG             0
#define WITH_CENTER_MASS    1
#define WITH_DARWIN         1

struct delta
{
    double rx,ry;
    double px,py;
};

struct particle
{
    double rx,ry;       // Position
    double px,py;       // Momentum
    float m;            // Mass
    float q;            // Charge
};

struct env
{
    int N;              // Number of particles
    struct particle *p; // Particle list
    struct particle *np;// Particle list with updated pos/momentum
    struct delta    *d;

    float M;            // Mass of the central object
    float dt;           // Time step
    float Rmin,         // Minimum and maximum disk radii
          Rmax;

    float t;            // Current simulation time
    float T;            // Total simulation time
    float Etot;         // Total energy

    float Coulomb_eps2; // Softening for the Coulomb potential

    float c;            // Speed of light
};

//------------------------------------------------------------------------------
// Macros
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Prototypes
//------------------------------------------------------------------------------

#endif 

