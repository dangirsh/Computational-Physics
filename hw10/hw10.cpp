/*
AEP 4380 HW #10
Dan Girshovich
4/25/13
The FFT and Spectral Methods
Compile with: g++ -std=c++0x -O2 -o gen_data hw10.cpp
Tested on Mac OSX 10.8.2 with a Intel Core 2 Duo
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "_array.hpp"
#include <fftw3.h>

using namespace std;

const int N = 256;
const double L = 300; // cm
const double delta = L / N; // cm
const double v = 300; // cm/s
const double x_A = 0.45 * L; // cm
const double y_A = 0.5 * L; // cm
const double s_A = 10; // cm
const double pi = 4 * atan(1.0);

fftw_complex *psi;
fftw_complex *d;
fftw_complex *dcos;
fftw_plan f_plan;
fftw_plan b_plan;

// returns the value of psi(x, y, t=0) given x, y
double psi_init(double x, double y) {
    double norm = (x - x_A) * (x - x_A) + (y - y_A) * (y - y_A);
    double gaussian = exp((-1) * norm / (s_A * s_A));
    // add square pulse
    if (x > (0.6 * L) && x < (0.7 * L) && y > (0.4 * L) && y < (0.5 * L)) {
        return gaussian + 1;
    } else {
        return gaussian;
    }
}

// sets psi to the inital value of gaussian + square pulse
void init_psi() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            psi[i * N + j][0] = psi_init(i * delta, j * delta);
        }
    }
}

// initalizes fftw data for computing ffts
void init_fft() {
    psi = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N * N);
    d = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N * N);
    dcos = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N * N);
    f_plan = fftw_plan_dft_2d(N, N, psi, d, FFTW_FORWARD, FFTW_ESTIMATE);
    b_plan = fftw_plan_dft_2d(N, N, dcos, psi, FFTW_BACKWARD, FFTW_ESTIMATE);
}

// returns |k_ij|
double mag_k(int i, int j) {
    if(i > (N / 2)) i = i - N;
    if(j > (N / 2)) j = j - N;
    double k_x = i / (N * delta), k_y = j / (N * delta);
    return 2 * pi * sqrt(k_x * k_x + k_y * k_y);
}

// returns d_ij * cos(v * t * |k_ij|)
void calc_dcos(double t) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            dcos[i * N + j][0] = d[i * N + j][0] * cos(v * t * mag_k(i, j));
            dcos[i * N + j][1] = d[i * N + j][1] * cos(v * t * mag_k(i, j));
        }
    }
}

// propagates the values in psi up to time t
void propagate_to(double t) {
    calc_dcos(t);
    fftw_execute(b_plan);
    // renormalize
    for(int i = 0; i < N * N; i++) psi[i][0] = psi[i][0] / (N * N);
}

// outputs the data file for psi at t
void write_psi(double t) {
    stringstream fname;
    fname << "dat/psi_" << t << ".dat";
    ofstream of;
    of.open(fname.str().c_str());
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            of << psi[j * N + i][0] << " ";
        }
        of << endl;
    }
    of.close();
}

int main() {
    init_fft();
    init_psi();
    // calculates and sets the d_ij coefficients for psi(x, t=0)
    fftw_execute(f_plan);
    double times[] = {0, 0.1, 0.2, 0.3};
    for (double t : times) {
        propagate_to(t);
        write_psi(t);
    }
}