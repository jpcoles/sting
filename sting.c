#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>

#include "sting.h"
#include "ic.h"
#include "log.h"
#include "fg.h"
#include "darwin.h"


int echolog=1;
FILE *logfp;

void help()
{
    fprintf(stderr, "Usage: sting\n");
    exit(2);
}

void energy(struct env *env, struct energy *e)
{
    int i,j;

    //double E = 0.0;
    double C, D=0;
    double Ceps2 = env->Coulomb_eps2;

    e->G =
    e->K =
    e->C = 
    e->D = 0;

    for (i=0; i < env->N; i++) 
    {
#if WITH_FG || WITH_CENTER_MASS
        e->G += -env->M / sqrt(pow(env->p[i].rx,2) + pow(env->p[i].ry,2) + env->grav_eps2);
#endif

        double p2 = pow(env->p[i].px,2) + pow(env->p[i].py,2);
        e->K += p2 / (2*env->p[i].m)
                * (1 - env->Kr*p2/(4*pow(env->p[i].m,2) * pow(env->c,2)));

        for (j=i+1; j < env->N; j++)
        {   
            //if (i==j) continue;

            double dx = env->p[i].rx - env->p[j].rx;
            double dy = env->p[i].ry - env->p[j].ry;

            double r2 = Ceps2 + dx*dx + dy*dy;
            double r  = sqrt(r2);

            C = env->p[i].q * env->p[j].q / r;

#if WITH_DARWIN
            D = env->p[i].px * env->p[j].px + env->p[i].py * env->p[j].py 
              + (env->p[i].px * dx + env->p[i].py * dy) * (env->p[j].px * dx + env->p[j].py * dy) / r2;

            //D /= 2 * env->p[i].m * env->p[j].m * pow(env->c,2);

//            + (env->p[i].px * dx - env->p[i].py * dy) * (env->p[j].px * dx - env->p[j].py * dy);
//          D /= r2 * (2 * env->p[i].m * env->p[j].m * pow(env->c,2));
#endif

            e->C += env->Kc * C;
            e->D -= env->Kd * env->Kc * C*D;
        }
    }

    e->E = e->G + e->K + e->C + e->D;
}

double ang_mom(struct env *env)
{
    int i;
    double J = 0;

    for (i=0; i < env->N; i++)
        J += env->p[i].rx * env->p[i].py - env->p[i].ry * env->p[i].px;

    return J;
}

double magnetic_field(struct env *env)
{
    int i;
    double J = 0;

    for (i=0; i < env->N; i++)
    {
        //fprintf(stderr, "%f %f\n", env->p[i].dHdpy, env->p[i].dHdpx);
        J += env->p[i].rx * env->p[i].dHdpy - env->p[i].ry * env->p[i].dHdpx;
    }

    return J;
}

inline void drift(struct env *env, double dt)
{
    int i;
#if WITH_FG
    for (i=0; i < env->N; i++) 
        fg(env->M, &env->p[i], dt);
#else
    for (i=0; i < env->N; i++) 
    {
        env->p[i].rx += dt * env->p[i].px / env->p[i].m;
        env->p[i].ry += dt * env->p[i].py / env->p[i].m;
    }
#endif
}

int main(int argc, char **argv)
{
    struct env env;
    struct energy E;

    set_units(&env);

    env.N    = 0;
    env.Rmin = .1;
    env.Rmax = 1.0;
    env.M    = 1e1;
    env.T    = 200;
    env.dt   = .01;
    env.Coulomb_eps2 = 0.05;
    env.grav_eps2 = 0.01;
    // env.c = 100;

    //--------------------------------------------------------------------------
    // Read command line arguments
    //--------------------------------------------------------------------------

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {"log", optional_argument, 0, 0},
        {"no-CM", no_argument, 0, 0},
        {0, 0, 0, 0}
    };

    if (argc < 2) help();

    while (1)
    {
        int option_index = 0;
        int c = getopt_long(argc, argv, "hvN:",
                            long_options, &option_index);

        if (c == -1) break;

        switch (c)
        {
            case 0:
                if (!strcmp("log", long_options[option_index].name))
                {
                    if (optarg)
                    {
                        if (open_log(optarg)) exit(1);
                    }
                }
                else if (!strcmp("no-CM", long_options[option_index].name))
                {
                    env.M = 0;
                }


                break;
            case 'h': help();      break;
            case 'N':
                env.N = atoi(optarg);
                if (env.N < 1)
                {
                    flog("Error: Number of particles must be more than 1.\n");
                    exit(1);
                }
                break;
            case '?': break;
            case ':': help(); break;
        }
    }

    if (logfp == NULL)
    {
        if (open_log("sting.log"))
        {
            fprintf(stderr, "Warning: Can't open default log file.\n");
        }
    }

    flog("sting version 0.1\n");

    //--------------------------------------------------------------------------
    // Create initial conditions
    //--------------------------------------------------------------------------

    srand48(0);

    fprintf(stderr, "sizeof(*p)  = %ld\n", sizeof(*env.p));
    fprintf(stderr, "sizeof(*np) = %ld\n", sizeof(*env.np));
    fprintf(stderr, "sizeof(*d)  = %ld\n", sizeof(*env.d));
    env.p  = malloc(env.N * sizeof(*env.p));  assert(env.p  != NULL);
    env.np = malloc(env.N * sizeof(*env.np)); assert(env.np != NULL);
    env.d  = malloc(env.N * sizeof(*env.d));  assert(env.d  != NULL);
    //ic_2_particle_simple(&env);
    //ic_4_particle_simple(&env);
    ic_random_circular(&env);
    //ic_random_elliptic(&env);
    //ic_coupled_pairs(&env);
    //ic_tube(&env);
    //ic_2tubes(&env);

    //--------------------------------------------------------------------------
    // During the interation loop, we will count the number of steps rather
    // than accumulate time in env.t in order to avoid round off errors for
    // large numbers of steps.
    //--------------------------------------------------------------------------
    int nsteps = ceil(env.T / env.dt);
    int step;
    int i,j;

    fprintf(stderr, "M = %f\n", env.M);
    fprintf(stderr, "eps2 = %f\n", env.Coulomb_eps2);
    for (step = 0; step < nsteps; step++)
    {
        env.t = step * env.dt;

        energy(&env, &E);
        env.Etot = E.E;
        env.Jtot = ang_mom(&env);
        env.Mtot = magnetic_field(&env);

        env.Mtot /= env.Jtot;

        for (i=0; i < env.N; i++)
        {
            printf("%i %i % 12e % 12e % 12e % 12e % 12e % 12e % 12f % 12e % 12e % 12e % 12e % 12e\n", 
                   step, i, 
                   env.p[i].rx, env.p[i].ry, env.p[i].px, env.p[i].py, env.p[i].q, env.Etot, env.Jtot, env.Mtot,
                   E.G, E.K, E.C, E.D);
        }
        printf("\n");
        //fprintf(stderr, "%f\n", env.Etot);

        drift(&env, env.dt/2);

        darwin_step(&env);
        //darwin_step(&env);

        drift(&env, env.dt/2);

    }



    return 0;
}

