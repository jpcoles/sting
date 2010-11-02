#include <stdlib.h>
#include <math.h>

#include "sting.h"
#include "ic.h"

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

        double px = 0;
        double py = v;

        env->p[i].m = 1.00;

        env->p[i].px = env->p[i].m * (px*cos(theta) - py*sin(theta));
        env->p[i].py = env->p[i].m * (px*sin(theta) + py*cos(theta));

        //-----------------------------------------------------------------------------
        // Assign a random charge (-1/+1)
        //-----------------------------------------------------------------------------
        env->p[i].q = (drand48() > 0.5) * 2 - 1;
    }

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
        double v = 0.5*sqrt(env->M / r);

        double px = 0;
        double py = v;

        env->p[i].m = 1.00;

        env->p[i].px = env->p[i].m * (px*cos(theta) - py*sin(theta));
        env->p[i].py = env->p[i].m * (px*sin(theta) + py*cos(theta));

        //-----------------------------------------------------------------------------
        // Assign a random charge (-1/+1)
        //-----------------------------------------------------------------------------
        env->p[i].q = (drand48() > 0.5) * 2 - 1;
    }

    return 0;
}
