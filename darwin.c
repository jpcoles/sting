#include <string.h>
#include "sting.h"
#include <math.h>

inline double dHdp(double mag2, double p, double m, double c)
{
    return -mag2 * p / (2 * pow(m,3) * pow(c,2));
}

#if 0
inline double dHdr(double mag2, double r, double m)
{
    return 
}
#endif

#if WITH_FG
int Hintegrator(struct env *env, struct particle *p, struct delta *d, double dt)
{
    int i,j;
    double pmag2;
    double dx, dy, r2, r;
    double dHdrx, dHdry;
    double darwpx, darwpy; 
    double dx0, dx1, dy0, dy1;
    double rhatx, rhaty;

    const double eps = 0.1;

    memset(d, 0, sizeof(*d) * env->N);

    for (i=0; i < env->N; i++)
    {
        pmag2 = pow(p[i].px,2) + pow(p[i].py,2);

        d[i].rx = 0;//dt * dHdp(pmag2, p[i].rx, p[i].m);
        d[i].ry = 0;//dt * dHdp(pmag2, p[i].ry, p[i].m);

        for (j=i+1; j < env->N; j++)
        {
            dx = p[i].rx - p[j].rx;
            dy = p[i].ry - p[j].ry;
            r2 = eps*eps + pow(dx, 2) + pow(dy, 2);
            r  = sqrt(r2);

            rhatx = dx/r;
            rhaty = dy/r;

           // d[i].rx += -p[i].q*p[j].q * (p[i].px+p[i].px*dx*dx/r2+dx*dy*p[i].py/r2) / (2*p[i].m*p[j].m*c*c*r);
           // d[i].ry += -p[i].q*p[j].q * (p[i].py+p[i].py*dy*dy/r2+dx*dy*p[i].px/r2) / (2*p[i].m*p[j].m*c*c*r);

            dx0 = p[i].py * p[j].px * rhaty
                + p[i].px * p[j].px * rhatx
                + p[i].px * p[j].py * rhaty
                - p[i].py * p[j].py * rhatx;

            dx1 = 3 * (  pow(rhatx,3) * p[i].px * p[j].px 
                       + pow(rhatx,2) * p[i].px * p[j].py * rhatx
                       + pow(rhatx,2) * p[i].py * p[j].px * rhaty
                       + pow(rhaty,2) * p[i].py * p[j].py * rhatx);

            dy0 = p[i].py * p[j].px * rhatx
                + p[i].px * p[j].py * rhatx
                - p[i].px * p[j].px * rhaty
                + p[i].py * p[j].py * rhaty;

            dy1 = 3 * (  pow(rhatx,2) * p[i].px * p[j].px * rhatx 
                       + pow(rhaty,2) * p[i].px * p[j].py * rhatx
                       + pow(rhaty,2) * p[i].py * p[j].px * rhatx
                       + pow(rhaty,3) * p[i].py * p[j].py);

            darwpx = p[i].q * p[j].q * (dx0-dx1) / (2 * p[i].m * p[j].m * c*c * r2);
            darwpy = p[i].q * p[j].q * (dy0-dy1) / (2 * p[i].m * p[j].m * c*c * r2);

            dHdrx = p[i].q * p[j].q * rhatx / r2; // + darwpx;
            dHdry = p[i].q * p[j].q * rhaty / r2; // + darwpy;

            //fprintf(stderr, "HI0 [%i,%i] %f %f  (%f %f) (%f %f) %f\n", i, j, dHdrx, dHdry, d[i].px, d[i].py, d[j].px, d[j].py, env->dt);

            d[i].px += dt * dHdrx;
            d[i].py += dt * dHdry;

            d[j].px -= dt * dHdrx;
            d[j].py -= dt * dHdry;

            //fprintf(stderr, "HI1 [%i,%i] %f %f  (%f %f) (%f %f) %f\n", i, j, dHdrx, dHdry, d[i].px, d[i].py, d[j].px, d[j].py, env->dt);
            //fprintf(stderr, "HI2 [%i,%i] %f %f  (%f %f) (%f %f) %f\n", i, j, dHdrx, dHdry, p[i].px, p[i].py, p[j].px, p[j].py, env->dt);
            //fprintf(stderr, "*-* %f %f %f %f %f %f %f %f %f %f %f \n", dx0, dx1, dy0, dy1, r2, dx, dy, dx*dx, dy*dy, p[i].rx, p[j].rx);
        }

        //fprintf(stderr, "!! %i %f %f %f %f \n", i, p[i].rx, p[i].ry, p[i].px, p[i].py);
    }

    return 0;
}

#else

int Hintegrator(struct env *env, struct particle *p, struct delta *d)
{
    int i,j;
    double pmag2, rmag2;
    double dx, dy, r2, r;
    double dHdrx, dHdry;
    double dHdpx, dHdpy;
    double darwpx, darwpy; 
    double dx0, dx1, dy0, dy1;
    double rhatx, rhaty;
    double R2, R;

    const double eps2 = env->Coulomb_eps2; //env->Rmax * 0.1;
    const double c = env->c;

    memset(d, 0, sizeof(*d) * env->N);

    for (i=0; i < env->N; i++)
    {
        pmag2 = pow(p[i].px,2) + pow(p[i].py,2);

        d[i].rx += dHdp(env->Kr*pmag2, p[i].rx, p[i].m, c);
        d[i].ry += dHdp(env->Kr*pmag2, p[i].ry, p[i].m, c);

#if WITH_CENTER_MASS
        R2 = env->grav_eps2 + pow(p[i].rx,2) + pow(p[i].ry,2);
        R  = sqrt(R2);
        d[i].px += -env->M*p[i].m * p[i].rx / (R2 * R);
        d[i].py += -env->M*p[i].m * p[i].ry / (R2 * R);
#endif

        p[i].dHdpx = 0;
        p[i].dHdpy = 0;

        for (j=0; j < env->N; j++)
        {
            if (i==j) continue;

            dx = p[i].rx - p[j].rx;
            dy = p[i].ry - p[j].ry;
            r2 = eps2 + pow(dx, 2) + pow(dy, 2);
            r  = sqrt(r2);

            rhatx = dx/r;
            rhaty = dy/r;

#if WITH_DARWIN
            dHdpx = (p[j].px + p[j].px*rhatx*rhatx + p[j].py*rhatx*rhaty);
            dHdpy = (p[j].py + p[j].px*rhatx*rhaty + p[j].py*rhaty*rhaty);

            //dHdpx /= 2*p[i].m*p[j].m*c*c;
            //dHdpy /= 2*p[i].m*p[j].m*c*c;

            dHdpx = env->Kd * env->Kc * p[i].q*p[j].q * dHdpx / r;
            dHdpy = env->Kd * env->Kc * p[i].q*p[j].q * dHdpy / r;

            p[i].dHdpx -= dHdpx; /* Doesn't include relativistic correction */
            p[i].dHdpy -= dHdpy;

            d[i].rx -= dHdpx;
            d[i].ry -= dHdpy;
#endif

            dx0 = p[i].px * p[j].px * rhatx + p[i].px * p[j].py * rhaty
                + p[i].px * p[j].px * rhatx + p[i].py * p[j].px * rhaty;

            dx1 = p[i].px * p[j].px * rhatx + p[i].py * p[j].py * rhatx
                + 3 * (  pow(rhatx,3) * p[i].px * p[j].px 
                       + pow(rhatx,2) * p[i].px * p[j].py * rhaty
                       + pow(rhatx,2) * p[i].py * p[j].px * rhaty
                       + pow(rhaty,2) * p[i].py * p[j].py * rhatx);

            dy0 = p[i].py * p[j].px * rhatx + p[i].py * p[j].py * rhaty
                + p[i].px * p[j].py * rhatx + p[i].py * p[j].py * rhaty;

            dy1 = p[i].px * p[j].px * rhaty + p[i].py * p[j].py * rhaty
                + 3 * (  pow(rhatx,2) * p[i].px * p[j].px * rhaty 
                       + pow(rhaty,2) * p[i].px * p[j].py * rhatx
                       + pow(rhaty,2) * p[i].py * p[j].px * rhatx
                       + pow(rhaty,3) * p[i].py * p[j].py);

            //darwpx = ; // / (2 * p[i].m * p[j].m * c*c);
            //darwpy = ; // / (2 * p[i].m * p[j].m * c*c);

            dHdrx = rhatx + env->Kd*(dx0-dx1);
            dHdry = rhaty + env->Kd*(dy0-dy1);

            //fprintf(stderr, "HI0 [%i,%i] %f %f  (%f %f) (%f %f) %f\n", i, j, dHdrx, dHdry, d[i].px, d[i].py, d[j].px, d[j].py, env->dt);

            d[i].px += env->Kc * p[i].q * p[j].q * dHdrx / r2;
            d[i].py += env->Kc * p[i].q * p[j].q * dHdry / r2;

            //d[j].px -= dHdrx;
            //d[j].py -= dHdry;

            //fprintf(stderr, "HI1 [%i,%i] %f %f  (%f %f) (%f %f) %f\n", i, j, dHdrx, dHdry, d[i].px, d[i].py, d[j].px, d[j].py, env->dt);
            //fprintf(stderr, "HI2 [%i,%i] %f %f  (%f %f) (%f %f) %f\n", i, j, dHdrx, dHdry, p[i].px, p[i].py, p[j].px, p[j].py, env->dt);
            //fprintf(stderr, "*-* %f %f %f %f %f %f %f %f %f %f %f \n", dx0, dx1, dy0, dy1, r2, dx, dy, dx*dx, dy*dy, p[i].rx, p[j].rx);
        }

        //fprintf(stderr, "!! %i %f %f \n", i, p[i].dHdpx, p[i].dHdpy);
        //fprintf(stderr, "!! %i %f %f %f %f \n", i, p[i].rx, p[i].ry, p[i].px, p[i].py);
    }


    return 0;
}

#endif

int darwin_step(struct env *env)
{
    int i, j, iter;

    struct particle *p  = env->p;
    struct particle *np = env->np;
    struct delta    *d  = env->d;

    memcpy(np, p, sizeof(*np) * env->N);

//  for (i=0; i < env->N; i++)
//      fprintf(stderr, "--- %i %f %f %f %f \n", i, p[i].rx, p[i].ry, p[i].px, p[i].py); 

    //--------------------------------------------------------------------------

    Hintegrator(env, np, d);

    //fprintf(stderr, "% 20.15f % 20.15f  % 20.15f  % 20.15f \n", d[0].rx, d[0].ry, d[0].px, d[0].py);

#if 1
    for (iter=0; iter < 10; iter++)
    {
        for (i=0; i < env->N; i++)
        {
            np[i].rx = p[i].rx + 0.5 * env->dt * d[i].rx;
            np[i].ry = p[i].ry + 0.5 * env->dt * d[i].ry;
            np[i].px = p[i].px + 0.5 * env->dt * d[i].px;
            np[i].py = p[i].py + 0.5 * env->dt * d[i].py;

        }

        Hintegrator(env, np, d);

//      for (i=0; i < env->N; i++)
//          fprintf(stderr, "## %i %f %f \n", i, p[i].dHdpx, p[i].dHdpy);

    }
#else
    for (i=0; i < env->N; i++)
    {
        d[i].px *= env->dt;
        d[i].py *= env->dt;
    }

#endif

    for (i=0; i < env->N; i++)
    {
            //fprintf(stderr, "## %i %f %f \n", i, p[i].dHdpx, p[i].dHdpy);
        p[i].rx = np[i].rx + 0.5 * env->dt * d[i].rx;
        p[i].ry = np[i].ry + 0.5 * env->dt * d[i].ry;
        p[i].px = np[i].px + 0.5 * env->dt * d[i].px;
        p[i].py = np[i].py + 0.5 * env->dt * d[i].py;
        p[i].dHdpx = np[i].dHdpx;
        p[i].dHdpy = np[i].dHdpy;

        //fprintf(stderr, "+++ %i %f %f %f %f \n", i, p[i].rx, p[i].ry, p[i].px, p[i].py); 
    }


    return 0;
}
