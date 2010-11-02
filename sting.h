#ifndef STING_H
#define STING_H

#include <stdio.h>

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

    float M;            // Mass of the central object
    float dt;           // Time step
    float Rmin,         // Minimum and maximum disk radii
          Rmax;

    float t;            // Current simulation time
    float T;            // Total simulation time
};

//------------------------------------------------------------------------------
// Macros
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Prototypes
//------------------------------------------------------------------------------

#endif 

