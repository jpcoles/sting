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

int main(int argc, char **argv)
{
    struct env env;

    env.N = 0;
    env.Rmin = 0.1;
    env.Rmax = 10.0;
    env.M = 100;
    env.T = 20;
    env.dt = 0.01;

    //--------------------------------------------------------------------------
    // Read command line arguments
    //--------------------------------------------------------------------------

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {"log", optional_argument, 0, 0},
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

    env.p  = malloc(env.N * sizeof(*env.p));  assert(env.p  != NULL);
    env.np = malloc(env.N * sizeof(*env.np)); assert(env.np != NULL);
    env.d  = malloc(env.N * sizeof(*env.d));  assert(env.d  != NULL);
    //ic_2_particle_simple(&env);
    ic_4_particle_simple(&env);
    //ic_random_circular(&env);
    //ic_random_elliptic(&env);

    //--------------------------------------------------------------------------
    // During the interation loop, we will count the number of steps rather
    // than accumulate time in env.t in order to avoid round off errors for
    // large numbers of steps.
    //--------------------------------------------------------------------------
    int nsteps = ceil(env.T / env.dt);
    int step;
    int i;

    for (step = 0; step < nsteps; step++)
    {
        env.t = step * env.dt;

        //for (i=0; i < env.N; i++) fg(env.M, &env.p[i], env.dt/2);

        darwin_step(&env);

        //for (i=0; i < env.N; i++) fg(env.M, &env.p[i], env.dt/2);

        for (i=0; i < env.N; i++)
            printf("%i] %f %f %f %f %f\n", i, env.p[i].rx, env.p[i].ry, env.p[i].px, env.p[i].py, env.p[i].q);
        printf("\n");


    }



    return 0;
}

