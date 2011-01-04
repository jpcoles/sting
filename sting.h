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

    double dHdpx, dHdpy;
};

struct energy
{
    double E;
    double G, K, C, D;
};

struct units
{
    double G,M,L,T,Q;
    double Myr, kpc, Msun;
    double year, AU;
};

struct env
{
    int N;              // Number of particles
    struct particle *p; // Particle list
    struct particle *np;// Particle list with updated pos/momentum
    struct delta    *d;

    float M;            // Mass of the central object
    float m;            // Mass of the particles
    float dt;           // Time step
    float Rmin,         // Minimum and maximum disk radii
          Rmax;

    float t;            // Current simulation time
    float T;            // Total simulation time
    float Etot;         // Total energy
    float Jtot;         // Total angular momentum
    float Mtot;         // Total magnetic field

    float Coulomb_eps2; // Softening for the Coulomb potential
    float grav_eps2;    // Softening for the gravitational potential

    float c;            // Speed of light

    float Kd, Kc;       // Coefficients for the darwin and coulomb parts
                        // of the hamiltonian.
    float Kr;
                                
    struct units units;
};

//------------------------------------------------------------------------------
// Macros
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Prototypes
//------------------------------------------------------------------------------

#endif 

