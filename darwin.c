
#include <string.h>
#include "sting.h"
#include <math.h>

const double c = 1000;

inline double dHdp(double mag2, double p, double m)
{
    return -mag2 * p / (2 * pow(m,3) * pow(c,2));
}

#if 0
inline double dHdr(double mag2, double r, double m)
{
    return 
}
#endif

int Hintegrator(struct env *env, struct particle *p, struct delta *d)
{
    int i,j;
    double pmag2;
    double dx, dy, r2, r;
    double dHdrx, dHdry;
    double darwpx, darwpy, dx0, dx1, dy0, dy1;

    memset(d, 0, sizeof(*d) * env->N);

    for (i=0; i < env->N; i++)
    {
        pmag2 = pow(p[i].px,2) + pow(p[i].py,2);

        //d[i].rx = 0.5 * env->dt * p[i].px / p[i].m;
        //d[i].ry = 0.5 * env->dt * p[i].py / p[i].m;

        d[i].rx = 0;//.5 * env->dt * dHdp(pmag2, p[i].rx, p[i].m);
        d[i].ry = 0;//.5 * env->dt * dHdp(pmag2, p[i].ry, p[i].m);

        for (j=i+1; j < env->N; j++)
        {
            dx = p[i].rx - p[j].rx;
            dy = p[i].ry - p[j].ry;
            r2 = pow(dx, 2) + pow(dy, 2);
            r  = sqrt(r2);

           // d[i].rx += -p[i].q*p[j].q * (p[i].px+p[i].px*dx*dx/r2+dx*dy*p[i].py/r2) / (2*p[i].m*p[j].m*c*c*r);
           // d[i].ry += -p[i].q*p[j].q * (p[i].py+p[i].py*dy*dy/r2+dx*dy*p[i].px/r2) / (2*p[i].m*p[j].m*c*c*r);

            dx0 = p[i].py*p[j].px*dy/r
                + p[i].px*p[j].px*dx/r
                + p[i].px*p[j].py*dy/r
                - p[i].py*p[j].py*dx/r;

            dx1 = 3*(pow(dx/r,3)*p[i].px*p[j].px 
                +    pow(dx/r,2)*p[i].px*p[j].py*dy/r
                +    pow(dx/r,2)*p[i].py*p[j].px*dy/r
                +    pow(dy/r,2)*p[i].py*p[j].py*dx/r);


            dy0 = p[i].py*p[j].px*dx/r
                + p[i].px*p[j].py*dx/r
                - p[i].px*p[j].px*dy/r
                + p[i].py*p[j].py*dy/r;

            dy1 = 3*(pow(dx/r,2)*p[i].px*p[j].px*dx/r 
                +    pow(dy/r,2)*p[i].px*p[j].py*dx/r
                +    pow(dy/r,2)*p[i].py*p[j].px*dx/r
                +    pow(dy/r,3)*p[i].py*p[j].py);

            darwpx = p[i].q*p[j].q*(dx0-dx1)/(2*p[i].m*p[j].m*c*c*r2);
            darwpy = p[i].q*p[j].q*(dy0-dy1)/(2*p[i].m*p[j].m*c*c*r2);

            dHdrx = p[i].q * p[j].q * dx / (r2*r);//+darwpx;
            dHdry = p[i].q * p[j].q * dy / (r2*r);//+darwpy;

            d[i].px += 0.5 * env->dt * dHdrx;
            d[i].py += 0.5 * env->dt * dHdry;

            d[j].px -= 0.5 * env->dt * dHdrx;
            d[j].py -= 0.5 * env->dt * dHdry;

            fprintf(stderr, "*-* %f %f %f %f %f %f %f %f %f %f %f \n", dx0, dx1, dy0, dy1, r2, dx, dy, dx*dx, dy*dy, p[i].rx, p[j].rx);
        }

        fprintf(stderr, "!! %i %f %f %f %f \n", i, p[i].rx, p[i].ry, p[i].px, p[i].py);
    }
    return 0;
}

int darwin_step(struct env *env)
{
    int i, j, iter;

    struct particle *p  = env->p;
    struct particle *np = env->np;
    struct delta    *d  = env->d;

    memcpy(np, p, sizeof(*np) * env->N);

    //--------------------------------------------------------------------------

    Hintegrator(env, p, d);

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

        Hintegrator(env, np, d);

        //fprintf(stderr, "% 20.15f % 20.15f  % 20.15f  % 20.15f \n", d[0].rx, d[0].ry, d[0].px, d[0].py);
    }
#endif

    env->Etot = 0.0;
    for (i=0; i < env->N; i++)
    {
        p[i].rx = np[i].rx + 0.5 * d[i].rx;
        p[i].ry = np[i].ry + 0.5 * d[i].ry;
        p[i].px = np[i].px + 0.5 * d[i].px;
        p[i].py = np[i].py + 0.5 * d[i].py;
        fprintf(stderr, "+++ %i %f %f %f %f \n", i, p[i].rx, p[i].ry, p[i].px, p[i].py); 
        env->Etot += - env->M*p[i].m/(hypot(p[i].rx,p[i].ry)) + (pow(p[i].px,2) + pow(p[i].py,2))/(2.0*p[i].m); 

        for (j=i+1; j < env->N; j++)
        {
            env->Etot += p[i].q * p[j].q / hypot(p[i].rx-p[j].rx, p[i].ry - p[j].ry);
        }
    }


    return 0;
}
