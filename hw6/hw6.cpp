/*
AEP 4380 HW #6
Dan Girshovich
3/28/13
Finds the potential and field inside a cylindrically symmetric electrode setup.
Compile with: g++ -std=c++0x -O3 -o gen_data hw6.cpp
Tested on Mac OSX 10.8.2 with a Intel Core 2 Duo
*/

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "limits.h"
#include <cmath>
#include <functional>
// #define ARRAYT_BOUNDS_CHECK
#include "array.hpp"

bool feq(float, float);
float part1(float, float);
void part2(float);
void part3(float, float, float);
void part4(float, float, float);
void part5(float, float, float);
float part6(float, float, float);

// units are mm
const float SPACE_TOL = 0.0001; // tol for counting grid pts as boundary pts
const float R_MIN = 0, R_MAX = 15, Z_MAX = 15, Z_MIN = -15;
const float Z1 = -2, Z_MID = 0, Z2 = 1.5, Z3 = 3, R1 = 2, R2 = 3, R3 = 10;
// units are V
const float V1 = 0, V2 = 3000, V3 = V1;

// floating point equality
bool feq(float a, float b) {
    return fabs(a - b) < SPACE_TOL;
};

class solution {
public:
    const int I_MAX, J_MAX;
    int num_its;
    solution(float _h, float omega, float dv_max);
    void for_each(std::function<void(int, int)> f);
    float i_to_z(int i);
    float j_to_r(int j);
    int z_to_i(float z);
    int r_to_j(float r);
    float get_v(float z, float r);
    float get_e_r(float z, float r);
    float get_e_z(float z, float r);
private:
    const float h; // mm
    array<bool> bound_flags; // true if point on boundary or electrode
    array<float> v; // values of the potential on the grid
    bool on_electrode_2(int i, int j);
    bool on_boundary(int i, int j);
    float laplace(int i, int j);
    float next_v(int i, int j, float omega);
    void solve(float omega, float dv_max);
};

// assumes valid spacing
solution::solution(float _h, float omega, float dv_max) :
    I_MAX((Z_MAX - Z_MIN) / _h), J_MAX((R_MAX - R_MIN) / _h),
    num_its(0), h(_h),
    bound_flags(I_MAX + 1, J_MAX + 1), v(I_MAX + 1, J_MAX + 1) {
    for_each([ = ](int i, int j) {
        bound_flags(i, j) = on_boundary(i, j);
        v(i, j) = on_electrode_2(i, j) ? V2 : V1;
    });
    solve(omega, dv_max);
}

// iterator for the arrays
void solution::for_each(std::function<void(int, int)> f) {
    for (int j = 0; j <= J_MAX; j++)
        for (int i = 0; i <= I_MAX; i++)
            f(i, j);
}

// --- Helpers ----
float solution::i_to_z(int i) {
    return Z_MIN + i * h;
}

float solution::j_to_r(int j) {
    return R_MIN + j * h;
}

int solution::z_to_i(float z) {
    return (z - Z_MIN) / h;
}

int solution::r_to_j(float r) {
    return (r - R_MIN) / h;
}

// --- Accessors ----
float solution::get_v(float z, float r) {
    return v(z_to_i(z), r_to_j(r));
}

// use numerical derivitive to find E_r(z, r) (V/m)
float solution::get_e_r(float z, float r) {
    int i = z_to_i(z), j = r_to_j(r);
    float e_r;
    if (j == 0) { // forward difference
        e_r = (v(i, j + 1) - v(i, j)) / h;
    } else if (j == J_MAX) { // backward difference
        e_r = (v(i, j) - v(i, j - 1)) / h;
    } else { // central difference
        e_r = (v(i, j + 1) - v(i, j - 1)) / (2 * h);
    }
    return -1000 * e_r; // V/m
}

// use numerical derivitive to find E_z(z, r) (V/m)
float solution::get_e_z(float z, float r) {
    int i = z_to_i(z), j = r_to_j(r);
    float e_z;
    if (i == 0) { // forward difference
        e_z = (v(i + 1, j) - v(i, j)) / h;
    } else if (i == I_MAX) { // backward difference
        e_z = (v(i, j) - v(i - 1, j)) / h;
    } else { // central difference
        e_z = (v(i + 1, j) - v(i - 1, j)) / (2 * h);
    }
    return -1000 * e_z; // V/m
}

// --- boundary helpers ---
bool solution::on_electrode_2(int i, int j) {
    float z = i_to_z(i), r = j_to_r(j);
    return feq(z, Z_MID) && r > R2 && r < R3;
}

bool solution::on_boundary(int i, int j) {
    float z = i_to_z(i), r = j_to_r(j);
    bool on_edge = feq(z, Z_MIN) || feq(z, Z_MAX) || feq(r, R_MAX);
    bool on_electrode_1 = feq(z, Z1) && r > R1 && r < R3;
    bool on_electrode_3 = z > Z2 && z < Z3 && r > R2 && r < R3;
    return on_edge || on_electrode_1 || on_electrode_2(i, j) || on_electrode_3;
}

// finite Laplace equation in cylindrical coordinates
float solution::laplace(int i, int j) {
    if (j == 0) { // special case on axis of symmetry
        return (1.0 / 6) * (4 * v(i, 1) + v(i + 1, 0) + v(i - 1, 0));
    } else {
        float a = v(i - 1, j), b = v(i, j + 1), c = v(i + 1, j), d = v(i, j - 1);
        return (a + b + c + d) / 4 + (b - d) / (8 * j);
    }
}

// performs an SOR step at i, j with omega
float solution::next_v(int i, int j, float omega) {
    float v_old = v(i, j);
    if (bound_flags(i, j)) { // electrodes & edges are at fixed V
        return v_old;
    } else {
        return omega * laplace(i, j) + (1.0 - omega) * v_old;
    }
}

// performs SOR with omega until the maximum change is < dv_max
void solution::solve(float omega, float dv_max) {
    float largest_dv = 0;
    do {
        num_its++;
        largest_dv = 0;
        for_each([ = , &largest_dv](int i, int j) {
            float v_new = next_v(i, j, omega);
            float dv = fabs(v(i, j) - v_new);
            if (dv > largest_dv) largest_dv = dv;
            v(i, j) = v_new;
        });
    } while (largest_dv > dv_max);
}

int main() {
    float best_omega = part1(0.1, 1);
    part2(best_omega);
    float h = 0.025, dv_max = 0.01; // from part2 data
    part3(h, best_omega, dv_max);
    part4(h, best_omega, dv_max);
    part5(h, best_omega, dv_max);
    printf("%f\n", part6(h, best_omega, dv_max));
}

// num_its vs. omega
float part1(float h, float dv_max) {
    return 1.85; // decided on this value
    float best_omega = 0.0;
    int smallest_num_its = INT_MAX;
    std::ofstream f;
    f.open("omega.dat");
    for (float omega = 1; omega < 1.99; omega += 0.05) {
        solution s(h, omega, dv_max);
        if (s.num_its < smallest_num_its) {
            smallest_num_its = s.num_its;
            best_omega = omega;
        }
        std::cout << omega << " " << s.num_its << std::endl;
    }
    f.close();
    return best_omega;
}

// finds v(r = 0, z = 0) for various h and dv_max
void part2(float omega) {
    std::ofstream f;
    f.open("accuracy_table.dat");
    for (float dv_max = 1; dv_max > 0.0001; dv_max /= 10) {
        for (float h = 0.1; h > 0.001; h /= 2) {
            solution s(h, omega, dv_max);
            f << h << " " << dv_max << " " << s.get_v(0, 0) << std::endl;
        }
    }
    f.close();
}

// Plots V and Ez as a function of z along r = 0
void part3(float h, float omega, float dv_max) {
    std::ofstream f_v, f_e_z;
    f_v.open("v_on_axis.dat");
    f_e_z.open("e_z_on_axis.dat");
    solution s(h, omega, dv_max);
    for (int i = 0; i <= s.I_MAX; i++) {
        float z = s.i_to_z(i);
        f_v << z << " " << s.get_v(z, 0) << std::endl;
        f_e_z << z << " " << s.get_e_z(z, 0) << std::endl;
    }
    f_v.close();
    f_e_z.close();
}

// Plots V  Er and Ez with r = 0.75mm. for −10mm. < z < +10mm
void part4(float h, float omega, float dv_max) {
    std::ofstream f_v, f_e_z, f_e_r;
    f_v.open("v_at_r75.dat");
    f_e_z.open("e_z_at_r75.dat");
    f_e_r.open("e_r_at_r75.dat");
    solution s(h, omega, dv_max);
    float r = 0.75;
    for (float z = -10; z < 10; z += h) {
        f_v << z << " " << s.get_v(z, r) << std::endl;
        f_e_z << z << " " << s.get_e_z(z, r) << std::endl;
        f_e_r << z << " " << s.get_e_r(z, r) << std::endl;
    }
    f_v.close();
    f_e_z.close();
    f_e_r.close();
}

// produces data for a contour plot of V near the electrodes and axis
void part5(float h, float omega, float dv_max) {
    std::ofstream f_v, f_z, f_r;
    f_v.open("v_contour.dat");
    f_z.open("z_contour.dat");
    f_r.open("r_contour.dat");
    solution s(h, omega, dv_max);
    // select range near electrodes and axis
    int lower_i = (s.I_MAX / 3), upper_i = 2 * (s.I_MAX / 3);
    int lower_j = 0, upper_j = s.J_MAX / 4;
    s.for_each([&](int i, int j) {
        if (i >= lower_i && i <= upper_i && j >= lower_j && j <= upper_j) {
            float z = s.i_to_z(i), r = s.j_to_r(j);
            f_v << s.get_v(z, r) << " ";
            if (i == upper_i) {
                f_r << r << std::endl;
                f_v << std::endl;
            }
            if (j == lower_j) f_z << z << std::endl;
        }
    });
    f_v.close();
    f_z.close();
    f_r.close();
}

// calculates the capacitance of the electrodes (pF)
float part6(float h, float omega, float dv_max) {
    float eps_0 = 8.8542, pi = 4.0 * atan(1.0), field_energy = 0.0;
    solution s(h, omega, dv_max);
    s.for_each([&](int i, int j) {
        float z = s.i_to_z(i), r = s.j_to_r(j);
        float e_z = s.get_e_z(z, r), e_r = s.get_e_r(z, r);
        field_energy += (e_z * e_z + e_r * e_r) * j;
    });
    field_energy *= pi * eps_0 * h * h * h * 1e-9; // mm -> m
    float v_cap = V2 - V1;
    return 2 * field_energy / (v_cap * v_cap);
}