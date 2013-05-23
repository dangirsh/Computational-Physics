/*
AEP 4380 HW #4
Dan Girshovich
2/28/13
Generates plot data and root info related to Bessel functions.
Compile with: g++ -std=c++0x -O2 -o gen_data hw3.cpp
Tested on Mac OSX 10.8.2 with a Intel Core 2 Duo
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "nr3.h"

// type / func defs
typedef void (*d_func)(const Doub, VecDoub &, VecDoub &);
void oscillator_derivs(const Doub, VecDoub &, VecDoub &);
void jc_derivs(const Doub, VecDoub &, VecDoub &);
void rk4(VecDoub &, VecDoub &, const Doub, const Doub, VecDoub &, d_func);
void gen_oscillator_data(const char *, double, double, int);
void gen_jc_data(const char *, double, double, int);
FILE *open_file_w(const char *);

const double w = 1.0; // oscillator frequency
const double mu = 1.0;
const double beta = 1.0;
const double pi = 4.0 * atan(1.0);

int main()
{
    gen_oscillator_data("oscillator_sol.dat", 0.0, 6 * pi, 1000);
    gen_jc_data("jc_sol_10000.dat", 0.0, 200.0, 10000);
}

// finds a solution to the simple harmonic oscillator equations between t1 and
// t2 using n steps of rk4.
void gen_oscillator_data(const char *fname, double t1, double t2, int n)
{
    FILE *fp = open_file_w(fname);
    const double h = (t2 - t1) / n;
    double t;
    VecDoub y(2);
    VecDoub dydt(2);
    // initial conditions
    y[0] = 1;
    y[1] = 0;
    oscillator_derivs(t1, y, dydt);
    // will hold values for next step
    VecDoub y_next(2);
    for (int i = 0; i < n; i++)
    {
        t = t1 + i * h;
        fprintf(fp, "%16.8f %16.8f %16.8f \n", t, y[0], dydt[0]);
        rk4(y, dydt, t, h, y_next, oscillator_derivs);
        //printf("%16.8f\n", y_next[0]);
        y = y_next; // uses overloaded assignment in NRvector
        oscillator_derivs(t, y, dydt);
    }
    fprintf(fp, "%16.8f %16.8f %16.8f \n", t2, y[0], dydt[0]);
    fclose(fp);
}

// returns the oscillator derivatives of each y at t in dydt
void oscillator_derivs(const Doub t, VecDoub &y, VecDoub &dydt)
{
    dydt[0] = y[1];
    dydt[1] = -1 * w * w * y[0];
}

// generates data for the plots of each parameter in the Janyes-Cummings model
void gen_jc_data(const char *fname, double t1, double t2, int n)
{
    FILE *fp = open_file_w(fname);
    const double h = (t2 - t1) / n;
    double t;
    VecDoub f(5); // x, y, z, e0, e1
    VecDoub dfdt(5); // x', y', z', e0', e1'
    // initial conditions
    f[0] = f[1] = 0; // x(0) = y(0) = 0
    f[2] = 1; // z(0) = 1
    f[3] = 1e-6; // e0(0) = 1e-6
    f[4] = 0; // e1(0) = e0'(0) = 0
    jc_derivs(t1, f, dfdt);
    // will hold values for next step
    VecDoub f_next(5);
    for (int i = 0; i < n; i++)
    {
        t = t1 + i * h;
        fprintf(fp, "%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f \n",
            t, f[0], f[1], f[2], f[3], dfdt[3]); // t, x, y, z, e0, e0'
        rk4(f, dfdt, t, h, f_next, jc_derivs);
        f = f_next; // uses overloaded assignment in NRvector
        jc_derivs(t1, f, dfdt);
    }
    // safer to print manually instead of looping past what the input specifies
    fprintf(fp, "%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f \n",
        t2, f[0], f[1], f[2], f[3], dfdt[3]); // t, x, y, z, e0, e0'
    fclose(fp);
}

// returns the Jaynes-Cummings derivatives of each y at t in dydt
void jc_derivs(const Doub t, VecDoub &f, VecDoub &dfdt)
{
    dfdt[0] = -1 * f[1]; // x' = - y
    dfdt[1] = f[0] + f[3] * f[2]; // y' = x + ez
    dfdt[2] = -1 * f[3] * f[1]; // z' = -ey
    dfdt[3] = f[4]; // e0' = e1
    dfdt[4] = beta * dfdt[1] - mu * mu * f[3]; // e1' = By' - u^2 e0
}

// Numerical Recipes p909 - slightly modified
void rk4(VecDoub &y, VecDoub &dydt, const Doub t, const Doub h,
         VecDoub &yout, d_func derivs)
{
    Int n = y.size();
    VecDoub dym(n), dyt(n), yt(n);
    Doub hh = h * 0.5;
    Doub h6 = h / 6.0;
    Doub th = t + hh;
    for (Int i = 0; i < n; i++) yt[i] = y[i] + hh * dydt[i];
    derivs(th, yt, dyt);
    for (Int i = 0; i < n; i++) yt[i] = y[i] + hh * dyt[i];
    derivs(th, yt, dym);
    for (Int i = 0; i < n; i++)
    {
        yt[i] = y[i] + h * dym[i];
        dym[i] += dyt[i];
    }
    derivs(t + h, yt, dyt);
    for (Int i = 0; i < n; i++)
        yout[i] = y[i] + h6 * (dydt[i] + dyt[i] + 2.0 * dym[i]);
}


// opens a new file for output
FILE *open_file_w(const char *fname)
{
    FILE *fp;
    fp = fopen(fname, "w+");
    if (NULL == fp)
    {
        printf("cannot open file\n");
        return (EXIT_SUCCESS);
    }
    return fp;
}