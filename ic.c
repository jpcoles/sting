#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "sting.h"
#include "ic.h"

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
        double v = sqrt(env->M / r);

        env->p[i].m = 1.0;

        double px = 0;
        double py = env->p[i].m*v;

        env->p[i].px = (px*cos(theta) - py*sin(theta));
        env->p[i].py = (px*sin(theta) + py*cos(theta));

        //-----------------------------------------------------------------------------
        // Assign a random charge (-1/+1)
        //-----------------------------------------------------------------------------
        double K = 10. / env->N;
        env->p[i].q = K; //*((drand48() > 0.5) * 2 - 1);
    }

    env->Coulomb_eps2 = pow(env->Rmax / sqrt(env->N) / 20, 2);

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
