#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "sting.h"
#include "ic.h"

int set_units(struct env *env)
{
    double year = env->units.year = 365 * 24 * 60 * 60;   // [s]
    double AU   = env->units.AU   = 1.4959787e11;         // [m]
    double Msun = env->units.Msun = 1.98892e30;           // [kg]

    double Gsi = 6.67300e-11;                             // [m^3 kg^-1 s^-2]
    double L   = env->units.L = 206265 * AU;                // [m]
    double T   = env->units.T = 100    * year;              // [s]
    double M   = env->units.M = 1e-6 * Msun;              // [kg]
    double Q   = env->units.Q = 10e10;                    // [C]
                 env->units.G = Gsi * (pow(L,-3) * M * pow(T,2));

    env->c = 3e10 * T / L;
    env->m  = 1. ; /// (sqrt(2) * env->c);

    env->m = 1;
    env->c = 1;
    //env->m  = 0.00001; //1. / (sqrt(2) * env->c);
    env->Kd = 1. / (2 * pow(env->m,2) * pow(env->c,2));
    //env->Kc = 1. / (4*M_PI*8.85e-12 * M * pow(L,3) / pow(T,2) / pow(Q,2));

    //env->Kd = env->Kd * 1e8;
    //env->Kc = 1;

    env->Kc  = 1e4 * 1e-4;
    env->Kd *= 1; //1e-4 * 1e8 * 1e-6;

    env->Kc = 1e-10;
    env->Kd = -1e-2 / env->Kc;

    //env->Kd = 1e-2;

    //env->Kd = 0;
//  env->Kc = 0;
    env->Kr = 0;


    fprintf(stderr, "unit G = %e\n", env->units.G);
    //fprintf(stderr, "env->M = %e\n", env->M);
    fprintf(stderr, "1 M = %e kg\n", env->M);
    fprintf(stderr, "1 T = %e s\n", T );
    fprintf(stderr, "1 L = %e m\n", L );
    fprintf(stderr, "1 Q = %e C\n", Q );
    fprintf(stderr, "m   = %e\n", env->m );
    fprintf(stderr, "c   = %e\n", env->c );
    fprintf(stderr, "Kd   = %e\n", env->Kd );
    fprintf(stderr, "Kc   = %e\n", env->Kc );
    //fprintf(stderr, "1 V = %e\n", (L/T) / (kpc/Myr));

    return 0;
}

int ic_2_particle_simple(struct env *env)
{
    assert(env->N == 2);

    env->p[0].rx = 1;
    env->p[0].ry = 0;
    env->p[0].px = 0;
    env->p[0].py = 0;
    env->p[0].m  = 1.0;
    env->p[0].q  = 1.0;

    env->p[1].rx = -1;
    env->p[1].ry = 0;
    env->p[1].px = 0;
    env->p[1].py = 0;
    env->p[1].m  = 1.0;
    env->p[1].q  = -1.0;

    return 0;
}

int ic_4_particle_simple(struct env *env)
{
    assert(env->N == 4);

    env->p[0].rx = 1;
    env->p[0].ry = 1;
    env->p[0].px = 0;
    env->p[0].py = 0;
    env->p[0].m  = 1.0;
    env->p[0].q  = -1.0;

    env->p[1].rx = -1;
    env->p[1].ry = 1;
    env->p[1].px = 0;
    env->p[1].py = 0;
    env->p[1].m  = 1.0;
    env->p[1].q  = -1.0;

    env->p[2].rx = 1;
    env->p[2].ry = -1;
    env->p[2].px = 0;
    env->p[2].py = 0;
    env->p[2].m  = 1.0;
    env->p[2].q  =  1.0;

    env->p[3].rx = -1;
    env->p[3].ry = -1;
    env->p[3].px = 0;
    env->p[3].py = 0;
    env->p[3].m  = 1.0;
    env->p[3].q  =  1.0;

    return 0;
}

int ic_random_circular(struct env *env)
{
    int i;

    assert((env->N & 1) == 0);

    int Q = -1;
    double K = 1.; // / env->N;

    for (i=0; i < env->N; i++)
    {
        //-----------------------------------------------------------------------------
        // Choose a random position in the disc
        //-----------------------------------------------------------------------------
        double r = (env->Rmax - env->Rmin)*sqrt(drand48()) + env->Rmin;
        double theta = 2*M_PI*drand48();

        env->p[i].rx = r * cos(theta);
        env->p[i].ry = r * sin(theta);

        //-----------------------------------------------------------------------------
        // Assume a circular orbit
        //-----------------------------------------------------------------------------
        double v = sqrt(r*r * env->M / pow(r*r + env->grav_eps2, 1.5));

        env->p[i].m = 1.; ///(env->c*sqrt(2.0));

        double px = 0;
        double py = env->p[i].m*v;

        env->p[i].px = (px*cos(theta) - py*sin(theta));
        env->p[i].py = (px*sin(theta) + py*cos(theta));

        //-----------------------------------------------------------------------------
        // Assign a random charge (-1/+1)
        //-----------------------------------------------------------------------------
        env->p[i].q = Q * K; //*((drand48() > 0.5) * 2 - 1);

        //Q *= -1;
    }

    //env->Coulomb_eps2 = pow(env->Rmax / sqrt(env->N), 2);

    return 0;
}

int ic_random_elliptic(struct env *env)
{
    int i;

    for (i=0; i < env->N; i++)
    {
        //-----------------------------------------------------------------------------
        // Choose a random position in the disc
        //-----------------------------------------------------------------------------
        double r = (env->Rmax - env->Rmin)*drand48() + env->Rmin;
        double theta = 2*M_PI*drand48();

        env->p[i].rx = r * cos(theta);
        env->p[i].ry = r * sin(theta);

        //-----------------------------------------------------------------------------
        // Assume a circular orbit
        //-----------------------------------------------------------------------------
        double v = (1.0 - drand48()) * sqrt(env->M / r);

        double px = 0;
        double py = v;

        env->p[i].m = 1.00;

        env->p[i].px = env->p[i].m * (px*cos(theta) - py*sin(theta));
        env->p[i].py = env->p[i].m * (px*sin(theta) + py*cos(theta));

        //-----------------------------------------------------------------------------
        // Assign a random charge (-1/+1)
        //-----------------------------------------------------------------------------
        env->p[i].q = 10.0*((drand48() > 0.5) * 2 - 1);
    }

    return 0;
}

int ic_coupled_pairs(struct env *env)
{
    assert(env->N == 2);

    double R = 1;
    double r = 0.1;
    double C = 0;
    double m = env->m;

    env->p[0].m  = m;
    env->p[0].q  = -1.0;
    env->p[0].rx = C + R + r;
    env->p[0].ry = 0;
    env->p[0].px = 0;
    double d = 2*r;
    env->p[0].py = sqrt(env->Kc * d*r/pow(d*d + env->Coulomb_eps2, 1.5));

    //env->p[0].py = m*sqrt(env->Kc/(4*env->p[0].m*sqrt(r*r + env->Coulomb_eps2)));
    //env->p[0].py = sqrt(1/(2*env->p[0].m*(2*r + sqrt(env->Coulomb_eps2))));

    env->p[1].m  = m;
    env->p[1].q  = 1.0;
    env->p[1].rx = C + R - r;
    env->p[1].ry = 0;
    env->p[1].px = 0;
    env->p[1].py = -env->p[0].py;

//  env->p[2].m  = m;
//  env->p[2].q  = -.2;
//  env->p[2].rx = C - R - r;
//  env->p[2].ry = 0;
//  env->p[2].px = 0;
//  env->p[2].py = 0; //0.2*sqrt(1/(2*env->p[2].m*0.2));

//  env->p[3].m  = m;
//  env->p[3].q  = .2;
//  env->p[3].rx = C - R + r;
//  env->p[3].ry = 0;
//  env->p[3].px = 0;
//  env->p[3].py = -env->p[2].py;

    return 0;
}

int ic_tube(struct env *env)
{
    int i=0;
    int Q = -1;
    double K = 1. / env->N;

    for (i=0; i < env->N; i++)
    {
        env->p[i].m  = 1.0;
        env->p[i].q  = 1;
        env->p[i].q = Q * K;
        env->p[i].rx = drand48() * 0.5 + 0.5;
        env->p[i].ry = drand48() * 0.5 - 0.25;
        env->p[i].px = -1;
        env->p[i].py = 0;

        Q *= -1;
    }


    return 0;
}

int ic_2tubes(struct env *env)
{
    assert((env->N & 1) == 0);

    int i=0;
    int Q = -1;
    double K = 1. / env->N;

    for (i=0; i < env->N/2; i++)
    {
        env->p[i].m  = 1.0;
        env->p[i].q  = 1;
        env->p[i].q = Q * K;
        env->p[i].rx = drand48() * 1.0 + 2.5;
        env->p[i].ry = drand48() * 0.25 - 0.12;
        env->p[i].px = -1;
        env->p[i].py = 0;

        Q *= -1;
    }

    for (; i < env->N; i++)
    {
        env->p[i].m  = 1.0;
        env->p[i].q  = 1;
        env->p[i].q = Q * K;
        env->p[i].rx = drand48() * 1.0 - 2.5;
        env->p[i].ry = drand48() * 0.25 + 0.12;
        env->p[i].px = 1;
        env->p[i].py = 0;

        Q *= -1;
    }


    return 0;
}

