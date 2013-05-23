// AEP 4380 Final Project
// Dan Girshovich
// 5/14/13
// Generates orbits and data about the new solutions to the 3 body problem
// Compile: clang++ -std=c++0x -stdlib=libc++ -g -o ../output/$1/main ../cpp/main.cpp
// Tested on Mac OSX 10.8.2 with a Intel Core 2 Duo


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include "nr3.h"
#include "stepper.h"
#include "odeint.h"
#include "stepperdopr853.h"
// #include "stepperdopr5.h"
#include <fftw3.h>
#include "ran.h"
#include <ctime>


string orbit_name;
int num_bodies;
double period;
VecDoub ics, masses;
Output out;

const double pi = 4.0 * atan(1.0);
const double G = 1;
const double tol = 1e-8; // atol and rtol
const int num_cycles = 5; // should be >= 1 for fft
// const bool dense = true;
const bool dense = false;
// ignored unless dense
const int points_per_cycle = 500;
const int num_points = num_cycles * points_per_cycle;

bool perturb = true;
double perturb_amount = 5e-5;
double perturb_prob = 0.2;
Ran rnd(time(NULL));

double mag(double a, double b) {
    return sqrt(a * a + b * b);
}

ofstream mass_of;

void perturb_masses(VecDoub &m) {
    for (int i = 0; i < num_bodies; i++) {
        double dm = 0.0;
        double r = rnd.doub();
        if(r < 0.3333){
            dm = -1 * perturb_amount;
        }
        else if(r < 0.6666){
            dm = perturb_amount;
        }
        if(rnd.doub() < perturb_prob){
            m[i] += dm;
        }
        if(m[i] < 1e-4){
            m[i] = 1e-4;
        }
        mass_of << m[i] << ' ';
    }
    mass_of << endl;
}

// functor for the 'derivs' calculation in the n body problem
struct rhs_grav {
    VecDoub m;
    double my_t = -1.0;
    rhs_grav(VecDoub mm, int nn) : m(mm) {}
    int x_index(int k) {
        return k;
    }
    int y_index(int k) {
        return k + num_bodies;
    }
    int vx_index(int k) {
        return k + 2 * num_bodies;
    }
    int vy_index(int k) {
        return k + 3 * num_bodies;
    }
    void set_vals(double &x, double &y, double &vx, double &vy, int k, VecDoub &w) {
        x = w[x_index(k)];
        y = w[y_index(k)];
        vx = w[vx_index(k)];
        vy = w[vy_index(k)];
    }
    // w should contain {x1,..., xn, y1,..., yn, vx1,..., vxn, vy1,..., vyn}
    void operator() (const double t, VecDoub &w, VecDoub &dwdt) {
        if(t > my_t){
            my_t = t;
            if (perturb) perturb_masses(m);
        }

        double xi, xj, yi, yj, vxi, vxj, vyi, vyj;

        for (int i = 0; i < num_bodies; i++) {
            set_vals(xi, yi, vxi, vyi, i, w);
            double sumx = 0, sumy = 0;
            for (int j = 0; j < num_bodies; j++) {
                if (i != j) {
                    set_vals(xj, yj, vxj, vyj, j, w);
                    double dx = xj - xi, dy = yj - yi;
                    double dist = mag(dx, dy);
                    double den = dist * dist * dist;
                    sumx += m[j] * dx / den;
                    sumy += m[j] * dy / den;
                }
            }
            dwdt[x_index(i)] = vxi;
            dwdt[y_index(i)] = vyi;
            dwdt[vx_index(i)] = G * sumx;
            dwdt[vy_index(i)] = G * sumy;
        }
    }
};

// uses odeint with rhs_grav to generate the positions of the three bodies
// results are saved in the global out object
void propagate(double t) {
    rhs_grav derivs(masses, num_bodies);
    double abtol = tol, rtol = tol, hinit = 1e-4, hmin = 0.0;
    Odeint<StepperDopr853<rhs_grav> > ode(ics, 0, t, abtol, rtol, hinit, hmin,
                                          out, derivs);
    ode.integrate();
}

// writes the orbit trajectories
void write_orbit_data() {
    ofstream of;
    of.open("../output/" + orbit_name + "/orbits.dat");
    for (int i = 0; i < out.count; i++) {
        of << out.xsave[i] << " ";
        for (int j = 0; j < 6; j++) {
            of << out.ysave[j][i] << " ";
        }
        of << endl;
    }
    of.close();
}

// returns the point on the shape space sphere for to the 3 body configuration
// return value is not scaled with the hyper-radius
// see http://suki.ipb.ac.rs/3body/info.php
string get_sphere_point(double q[]) {
    double x1 = q[0], x2 = q[1], x3 = q[2], y1 = q[3], y2 = q[4], y3 = q[5];
    double rho_x = (1 / sqrt(2.0)) * (x1 - x2);
    double rho_y = (1 / sqrt(2.0)) * (y1 - y2);
    double rho = mag(rho_x, rho_y);
    double lambda_x = (1 / sqrt(6.0)) * (x1 + x2 - 2 * x3);
    double lambda_y = (1 / sqrt(6.0)) * (y1 + y2 - 2 * y3);
    double lambda = mag(lambda_x, lambda_y);
    double n_x = 2 * (rho_x * lambda_x + rho_y * lambda_y);
    double n_y = lambda * lambda - rho * rho;
    double n_z = 2 * (rho_x * lambda_y - rho_y * lambda_x);
    stringstream ss;
    ss << n_x << " " << n_y << " " << n_z;
    return ss.str();
}

// writes the shape sphere data
void write_sphere_data() {
    ofstream of;
    of.open("../output/" + orbit_name + "/sphere.dat");
    double q[6];
    for (int i = 0; i < out.count; i++) {
        for (int j = 0; j < 6; j++) q[j] = out.ysave[j][i];
        of << get_sphere_point(q) << endl;
    }
    of.close();
}

// writes the positions and velocities after every period
void write_phase_data() {
    ofstream of;
    of.open("../output/" + orbit_name + "/phases.dat");
    for (int i = 0, j = 0; i < num_cycles; i++, j += points_per_cycle) {
        for (int k = 0; k < 6; k++) {
            of << out.ysave[k][j] << " " << out.ysave[k + 6][j] << " ";
        }
        of << endl;
    }
    of.close();
}

// updates the array so that the values are centered around 0.0
void center(fftw_complex *q, int len) {
    double sum = 0.0;
    for (int i = 0; i < num_points; i++) sum += q[i][0];
    double average = sum / num_points;
    for (int i = 0; i < num_points; i++) q[i][0] = q[i][0] - average;
}


void write_fft() {
    if (!dense) {
        cout << "Error: fft requires dense output" << endl;
        return;
    }
    fftw_complex *fft_in, *fft_out;
    fftw_plan plan;
    fft_in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * num_points);
    fft_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * num_points);
    plan = fftw_plan_dft_1d(num_points, fft_in, fft_out, FFTW_FORWARD,
                            FFTW_ESTIMATE);
    for (int i = 0; i < num_points; i++) {
        fft_in[i][0] = out.ysave[0][i];
        fft_in[i][1] = 0; // all real
    }

    center(fft_in, num_points);

    fftw_execute(plan);

    ofstream of;
    string fname = "../output/" + orbit_name + "/fft_x1" + ".dat";
    of.open(fname);
    for (int i = 0; i < num_points; i++) {
        of << fft_out[i][0] << endl;
    }
    of.close();
}

// parses command line args into global vars
void init_globals(char **argv) {
    orbit_name = string(argv[1]);
    num_bodies = atoi(argv[2]);
    period = atof(argv[3]);
    masses = VecDoub(num_bodies);
    for (int i = 0; i < num_bodies; i++) masses[i] = atof(argv[i + 4]);
    ics = VecDoub(2 * 2 * num_bodies); // 2 dimensions * 2 eq/dim * # masses
    for (int i = 0; i < 4 * num_bodies; i++) ics[i] = atof(argv[num_bodies + 4 + i]);
    out = Output(dense ? num_points : 0);
    mass_of.open("../output/" + orbit_name + "/mass.dat");
}

int main(int argc, char **argv) {
    // TODO: check valid args
    init_globals(argv);
    propagate(num_cycles * period);
    write_orbit_data();
    // write_sphere_data();
    // write_phase_data();
    // write_fft();
    mass_of.close();
}