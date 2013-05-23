/*
AEP 4380 HW #4
Dan Girshovich
4/4/13
Solves the time dependent Schrodinger Eq for the propagation of an electron wave
Compile with: g++ -std=c++0x -O2 -o gen_data hw7.cpp
Tested on Mac OSX 10.8.2 with a Intel Core 2 Duo
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <functional>
// #define ARRAYT_BOUNDS_CHECK
#include "array.hpp"

typedef std::complex<double> cx;

// general constants
const cx i = cx(0, 1); // imaginary constant i
const double h_bar = 6.5821e-16; // eV-sec (Planck's constant / 2 pi)
const double k = 3.801; // eV-A^2 (h_bar^2 / 2m for the electron)

// for the potential
const double L = 1000; // A (length of region)
const double V_0 = 4.05; // eV (height of potential barrier)
const double omega_v = 7; // A (softness of potential barrier)

// for the wavefunction
const double s = 20; // A (width of Gaussian wave function)
const double k_0 = 1; // A^-1 (average wavenumber)

// for the propagation
const double dx = 0.1; // A sample size in space
const double dt = 0.01; // fs sample size in time
const double omega = 2 * dx * dx / (dt * 1e-15);
const int j_max = L / dx;
array<double> V = array<double>(j_max); // V(x) potential
array<cx> psi = array<cx>(j_max); // psi(x) wavefunction at some time

// elements in matrix equation for propagation
array<cx> a = array<cx>(j_max);
array<cx> b = array<cx>(j_max);
array<cx> c = array<cx>(j_max);
array<cx> d = array<cx>(j_max);

// psi is held at 0 at the ends
cx get_psi(int j) {
    return (j == -1 || j == j_max) ? cx(0, 0) : psi(j);
}

// sets the values for d based on psi's current state (in time)
void set_d() {
    for (int j = 0; j < j_max; j++) {
        cx z = (omega * i * h_bar + 2 * k + dx * dx * V(j)) / k;
        d(j) = -1.0 * get_psi(j - 1) + z * get_psi(j) - get_psi(j + 1);
    }
}

// gets a, b, c, d, V, and psi ready for propagation
void init_arrays() {
    for (int j = 0; j < j_max; j++) {
        double x = j * dx;
        V(j) = V_0 / (1 + exp((0.5 * L - x) / omega_v));
        double z = (x - 0.3 * L) / s;
        psi(j) = exp(-1 * z * z + i * x * k_0); // exp is overloaded for complex
        a(j) = cx(1, 0);
        b(j) = (omega * i * h_bar - 2 * k - dx * dx * V(j)) / k;
        c(j) = cx(1, 0);
        set_d();
    }
}

// From section 2.4 of Numerical Recipes. Modified to use template for matrix
// types and bounds checking array class.
template<class T>
void tridag(array<T> &a, array<T> &b, array<T> &c, array<T> &d, array<T> &u) {
    int j, n = a.n();
    T bet;
    array<T> gam(n);
    if (b(0) == 0.0) throw("Error 1 in tridag");
    u(0) = d(0) / (bet = b(0));
    for (j = 1; j < n; j++) {
        gam(j) = c(j - 1) / bet;
        bet = b(j) - a(j) * gam(j);
        if (bet == 0.0) throw("Error 2 in tridag");
        u(j) = (d(j) - a(j) * u(j - 1)) / bet;
    }
    for (j = (n - 2); j >= 0; j--)
        u(j) -= gam(j + 1) * u(j + 1);
}

void output_potential() {
    std::ofstream f;
    f.open("potential.dat");
    for (int j = 0; j < j_max; j++)
        f << j *dx << " " << V(j) << "\n";
    f.close();
}

void output_psi(std::string fname) {
    std::ofstream f;
    f.open(fname.c_str());
    for (int j = 0; j < j_max; j++) {
        cx p = psi(j);
        double x = j * dx;
        f << x << " " << p.real() << " " << p.imag() << " " << norm(p) << "\n";
    }
    f.close();
}

// simple Reimann sum to find the total probability
double integral_norm_psi() {
    double integral = 0;
    for (int j = 0; j < j_max; j++)
        integral += norm(psi(j)) * dx;
    return integral;
}

int main() {
    init_arrays();
    output_potential();
    double t = 0;
    for (int n = 0; n <= 5; n++) {
        std::stringstream fname;
        fname << "psi" << n << ".dat";
        output_psi(fname.str());
        std::cout << integral_norm_psi() << "\n";
        // propagation
        for (double t_next = t + 10; t < t_next; t += dt) {
            tridag(a, b, c, d, psi);
            set_d();
        }
    }
}