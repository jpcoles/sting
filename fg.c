#include <math.h>
#include <stdio.h>
#include "sting.h"


double dot(double *x0, double *x1)
{
    return x0[0]*x1[0] + x0[1]*x1[1];
}

void fg(double mu, struct particle *p, double dt)
{
    double          f, g, fd, gd;
    double          r, vsq, u, a, ec, es, e, en, nf;
    double          dec, dm, s, c, w, wp, wpp, wppp, dx;
    int             iter;
    double          lo, up, next;

    //r = sqrt(dot(x, x));
    //vsq = dot(v, v);
    //u = dot(x, v);
    r = hypot(p->rx, p->ry);
    vsq = pow(p->px,2) + pow(p->py,2);
    u = p->rx * p->px + p->ry * p->py;

    a = 1 / (2 / r - vsq / mu);
    en = sqrt(mu / (a * a * a));
    ec = 1 - r / a;
    es = u / (en * a * a);
    e = sqrt(ec * ec + es * es);

    nf = en / (2 * M_PI);
    dt -= ((int) (nf * dt)) / nf;
    dm = en * dt - es;

    if ((es * cos(dm) + ec * sin(dm)) > 0)
        dec = dm + 0.85 * e;
    else
        dec = dm - 0.85 * e;

    lo = -2 * M_PI;
    up = 2 * M_PI;
    for (iter = 1; iter <= 32; iter++) 
    {
        s = sin(dec);
        c = cos(dec);
        w = dec - ec * s - es * c - dm;

        if (w > 0)
            up = dec;
        else
            lo = dec;

        wp = 1 - ec * c + es * s;
        wpp = ec * s + es * c;
        wppp = ec * c - es * s;
        dx = -w / wp;
        dx = -w / (wp + dx * wpp / 2);
        dx = -w / (wp + dx * wpp / 2 + dx * dx * wppp / 6);
        next = dec + dx;

        if (dec == next)
            break;

        if ((next > lo) && (next < up))
            dec = next;
        else
            dec = (lo + up) / 2;

        if ((dec == lo) || (dec == up))
            break;

        if (iter > 30)
            printf("%4d %23.20f %e\n", iter, dec, up - lo);
    }

    if (iter > 32)
        printf("dec soln failed %f\n", 42.);

    f = (a / r) * (c - 1) + 1;
    g = dt + (s - dec) / en;
    fd = -(a / (r * wp)) * en * s;
    gd = (c - 1) / wp + 1;

//  for (i = 1; i <= 3; i++) 
//  {
//      w = f * x[i] + g * v[i];
//      v[i] = fd * x[i] + gd * v[i];
//      x[i] = w;
//  }

    w = f * p->rx + g * p->px;
    p->px = fd * p->rx + gd * p->px;
    p->rx = w;

    w = f * p->ry + g * p->py;
    p->py = fd * p->ry + gd * p->py;
    p->ry = w;
}
