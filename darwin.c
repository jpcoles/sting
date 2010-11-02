#include <math.h>
#include <string.h>
#include "sting.h"

inline double dHdp(double mag2, double p, double m)
{
    const double c = 1000;
    return -mag2 * p / (2 * pow(m,3) * pow(c,2));
}

#if 0
inline double dHdr(double mag2, double r, double m)
{
    return 
}
#endif

int darwin_step(struct env *env)
{
    int i,j, iter;

    struct particle *p  = env->p;
    struct particle *np = env->np;
    struct delta    *d  = env->d;

    double pmag2;
    double dx, dy, r2, r;
    double dHdrx, dHdry;

    memset(d, 0, sizeof(*env->d) * env->N);

    //--------------------------------------------------------------------------

    for (i=0; i < env->N; i++)
    {
        pmag2 = pow(p[i].px,2) + pow(p[i].py,2);

        //d[i].rx = 0.5 * env->dt * p[i].px / p[i].m;
        //d[i].ry = 0.5 * env->dt * p[i].py / p[i].m;

        d[i].rx = 0.5 * env->dt * dHdp(pmag2, p[i].rx, p[i].m);
        d[i].ry = 0.5 * env->dt * dHdp(pmag2, p[i].ry, p[i].m);

        for (j=i+1; j < env->N; j++)
        {
            dx = p[i].rx - p[j].rx;
            dy = p[i].ry - p[j].ry;
            r2 = pow(dx, 2) + pow(dy, 2);
            r  = sqrt(r2);

            dHdrx = p[i].q * p[j].q * dx / (r2*r);
            dHdry = p[i].q * p[j].q * dy / (r2*r);

            d[i].px += 0.5 * env->dt * dHdrx;
            d[i].py += 0.5 * env->dt * dHdry;

            d[j].px -= 0.5 * env->dt * dHdrx;
            d[j].py -= 0.5 * env->dt * dHdry;
        }

        //fprintf(stderr, "!! %i %f %f\n", i, d[i].px, d[i].py);
    }

    //--------------------------------------------------------------------------

#if 1
    for (iter=0; iter < 10; iter++)
    {
        for (i=0; i < env->N; i++)
        {
            np[i].rx = p[i].rx + 0.5 * d[i].rx;
            np[i].ry = p[i].ry + 0.5 * d[i].ry;
            np[i].px = p[i].px + 0.5 * d[i].px;
            np[i].py = p[i].py + 0.5 * d[i].py;
        }

        memset(d, 0, sizeof(*env->d) * env->N);

        for (i=0; i < env->N; i++)
        {
            pmag2 = pow(np[i].px,2) + pow(np[i].py,2);

            //d[i].rx = 0.5 * env->dt * np[i].px / p[i].m;
            //d[i].ry = 0.5 * env->dt * np[i].py / p[i].m;

            d[i].rx = 0.5 * env->dt * dHdp(pmag2, np[i].rx, p[i].m);
            d[i].ry = 0.5 * env->dt * dHdp(pmag2, np[i].ry, p[i].m);

            for (j=i+1; j < env->N; j++)
            {
                dx = np[i].rx - np[j].rx;
                dy = np[i].ry - np[j].ry;
                r2 = pow(dx, 2) + pow(dy, 2);
                r  = sqrt(r2);

                dHdrx = p[i].q * p[j].q * dx / (r2*r);
                dHdry = p[i].q * p[j].q * dy / (r2*r);

                d[i].px += 0.5 * env->dt * dHdrx;
                d[i].py += 0.5 * env->dt * dHdry;

                d[j].px -= 0.5 * env->dt * dHdrx;
                d[j].py -= 0.5 * env->dt * dHdry;
            }
        }

        fprintf(stderr, "% 20.15f % 20.15f  % 20.15f  % 20.15f \n", d[0].rx, d[0].ry, d[0].px, d[0].py);
    }

#endif

    for (i=0; i < env->N; i++)
    {
        p[i].rx = np[i].rx + 0.5 * d[i].rx;
        p[i].ry = np[i].ry + 0.5 * d[i].ry;
        p[i].px = np[i].px + 0.5 * d[i].px;
        p[i].py = np[i].py + 0.5 * d[i].py;
    }

    return 0;
}
