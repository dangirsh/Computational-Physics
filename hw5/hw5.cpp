/*
AEP 4380 HW #5
Dan Girshovich
3/7/13
Calculates the positions of three bodies acting under gravity.
Compile with: g++ -std=c++0x -O2 -o gen_data hw5.cpp
Tested on Mac OSX 10.8.2 with a Intel Core 2 Duo
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "nr3.h"
#include "stepper.h"
#include "odeint.h"
#include "stepperdopr5.h"

// type / func defs
void gen_traj_data(Doub, const char *);
void alert_if_collision(Int, Int, Doub, Doub);
FILE *open_file_w(const char *);

const Doub G = 6.6726e-11; // Nm^2/kg^s (Gravitational constant)

const Doub Me = 5.976e24; // kg (Mass of the Earth)
const Doub Mm = 0.0123 * Me; // kg (Mass of the Moon)
const Doub Mm2 = 0.2 * Mm; // kg (Mass of the second Moon)

const Doub Re = 6.387e6; // m (Radius of the Earth)
const Doub Rm = 3.476e6; // m (Radius of the Earth)
const Doub Rm2 = 0.5 * Rm; // m (Radius of the Earth)

int main()
{
    gen_traj_data(200, "out.dat");
}

// functor for the 'derivs' calculation in the n body problem
struct rhs_grav
{
    Int n; // number of bodies
    VecDoub m; // for the n masses
    rhs_grav(VecDoub mm, Int nn) : m(mm), n(nn) {}
    Int x_index(Int k) { return k; }
    Int y_index(Int k) { return k + n; }
    Int vx_index(Int k) { return k + 2 * n; }
    Int vy_index(Int k) { return k + 3 * n; }
    void set_vals(Doub &x, Doub &y, Doub &vx, Doub &vy, Int k, VecDoub &w)
    {
        x = w[x_index(k)];
        y = w[y_index(k)];
        vx = w[vx_index(k)];
        vy = w[vy_index(k)];
    }
    // w should contain {x1,..., xn, y1,..., yn, vx1,..., vxn, vy1,..., vyn}
    void operator() (const Doub t, VecDoub &w, VecDoub &dwdt)
    {
        Doub xi, xj, yi, yj, vxi, vxj, vyi, vyj;
        for (Int i = 0; i < n; i++)
        {
            set_vals(xi, yi, vxi, vyi, i, w);
            Doub sumx = 0, sumy = 0;
            for (Int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    set_vals(xj, yj, vxj, vyj, j, w);
                    Doub d = pow(pow(xi - xj, 2) + pow(yi - yj, 2), 1.5);
                    alert_if_collision(i, j, d, t);
                    sumx += m[j] * (xj - xi) / d;
                    sumy += m[j] * (yj - yi) / d;
                }
            }
            dwdt[x_index(i)] = vxi;
            dwdt[y_index(i)] = vyi;
            dwdt[vx_index(i)] = G * sumx;
            dwdt[vy_index(i)] = G * sumy;
        }
    }
};

// sets the initial positions and velocities for the 3 bodies
void set_ics(VecDoub &ic)
{
    // initial x values
    ic[0] = 0.0;
    ic[1] = 0.0;
    ic[2] = -6.0e8;
    //initial y values
    ic[3] = 0.0;
    ic[4] = 3.84e8;
    ic[5] = 4.05e8;
    // initial vx values
    ic[6] = -12.593;
    ic[7] = 1019.0;
    ic[8] = 500.0;
    // initial vy values
    ic[9] = 0.0;
    ic[10] = 0.0;
    ic[11] = 50.0;
}

// uses odeint with rhs_grav to generate the positions of the three bodies
void gen_traj_data(Doub days, const char *fname)
{
    Doub secs = days * 24 * 60 * 60;
    Int nbodies = 3;
    Output out(0);
    VecDoub m(3);
    // masses of Earth, Moon, and Moon2
    m[0] = Me, m[1] = Mm, m[2] = Mm2;
    VecDoub ic(nbodies * 4); // there are 2 * 2 * N initial conditions
    set_ics(ic);
    rhs_grav derivs(m, nbodies);
    Doub rtol = 1e-8, abtol = 0.0, hinit = 1.0, hmin = 0.0;
    Odeint<StepperDopr5<rhs_grav> > ode(ic, 0, secs, abtol, rtol, hinit, hmin,
                                        out, derivs);
    ode.integrate();
    FILE *fp = open_file_w(fname);
    // print calculated positions of every body for each time step
    for (Int i = 0; i < out.count; i++)
    {
        fprintf(fp, "%f %f %f %f %f %f %f \n",
                out.xsave[i], out.ysave[0][i], out.ysave[1][i], out.ysave[2][i],
                out.ysave[3][i], out.ysave[4][i], out.ysave[5][i]);
    }
    fclose(fp);
}

// debugging tool for detecting collisions
void alert_if_collision(Int i, Int j, Doub d, Doub t)
{
    VecDoub r(3);
    r[0] = Re, r[1] = Rm, r[2] = Rm2;
    if (r[i] + r[j] > d)
    {
        printf("Collision: %f %d %d \n", t, i, j);
        exit(1);
    }
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