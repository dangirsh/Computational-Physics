/*
AEP 4380 HW #8
Dan Girshovich
4/12/13
Least Squares Curve Fitting
Compile with: g++ -std=c++0x -O2 -o gen_data hw8.cpp
Tested on Mac OSX 10.8.2 with a Intel Core 2 Duo
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <functional>
#include "nr3.h"

using namespace std;

typedef function<double(double)> math_func;
const int N = 607; // max number of data points
const double sigma_perc =  0.2;
vector<int> ts = vector<int>(N);
vector<double> ys = vector<double>(N);
vector<double> sigmas = vector<double>(N);

// slightly modified given code to read from the data file
void load_arrays() {
    const int NCMAX = 200; // max number of characters per line
    char cline[N];
    double co2;
    FILE *fp;
    fp = fopen("maunaloa.co2", "r");
    if (NULL == fp) {
        printf("Can't open file.");
        exit(0);
    }
    for (int i = 0; i < 15; i++) fgets(cline, NCMAX, fp); // read a whole line
    int t = 0, npts = 0;
    int year;
    do {
        fscanf(fp, "%d", &year);
        for (int j = 0; j < 12; j++) {
            fscanf(fp, "%lf", &co2);
            if (co2 > 0.0) {
                ts[npts] = t;
                ys[npts++] = co2;
            }
            t++;
        }
        if (npts >= N) break;
        fscanf(fp, "%lf", &co2); // skipping average value
    } while (year < 2008);
}

// returns an array of fit functions
vector<math_func> get_fit_funcs(bool first_harmonic) {
    const double pi = 4.0 * atan(1.0);
    double k = 2 * pi / 12;
    vector<math_func> fs {
        [ = ](double t){return 1;},
        [ = ](double t){return t;},
        [ = ](double t){return t * t;},
        [ = ](double t){return sin(k * t);},
        [ = ](double t){return cos(k * t);},
    };
    vector<math_func> half_year_terms {
        [ = ](double t) {return sin((k * t) / 2);},
        [ = ](double t) {return cos((k * t) / 2);},
    };
    if (first_harmonic) {
        fs.insert(fs.end(), half_year_terms.begin(), half_year_terms.end());
    }
    return fs;
}

// returns the F matrix given fit functions
MatDoub_IO get_F(vector<math_func> fs) {
    int m = fs.size();
    MatDoub_IO F(m, m);
    for (int l = 0; l < m; l++) {
        for (int k = 0; k < m; k++) {
            double sum = 0.0;
            for (int i = 0; i < N; i++) {
                double y = ys[i];
                sum += fs[l](ts[i]) * fs[k](ts[i]) / (sigmas[i] * sigmas[i]);
            }
            F[l][k] = sum;
        }
    }
    return F;
}

// returns the b vector given fit functions
MatDoub_IO get_b(vector<math_func> fs) {
    int m = fs.size();
    MatDoub_IO b(m, 1);
    for (int l = 0; l < m; l++) {
        double sum = 0.0;
        for (int i = 0; i < N; i++) {
            double y = ys[i];
            sum += y * fs[l](ts[i]) / (sigmas[i] * sigmas[i]);
        }
        b[l][0] = sum;
    }
    return b;
}

// From Numerical Recipes Ch 2
void gaussj(MatDoub_IO &a, MatDoub_IO &b) {
    Int i, icol, irow, j, k, l, ll, n = a.nrows(), m = b.ncols();
    Doub big, dum, pivinv;
    VecInt indxc(n), indxr(n), ipiv(n);
    for (j = 0; j < n; j++) ipiv[j] = 0;
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++)
            if (ipiv[j] != 1)
                for (k = 0; k < n; k++) {
                    if (ipiv[k] == 0) {
                        if (abs(a[j][k]) >= big) {
                            big = abs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l = 0; l < n; l++) SWAP(a[irow][l], a[icol][l]);
            for (l = 0; l < m; l++) SWAP(b[irow][l], b[icol][l]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) throw("gaussj: Singular Matrix");
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (l = 0; l < n; l++) a[icol][l] *= pivinv;
        for (l = 0; l < m; l++) b[icol][l] *= pivinv;
        for (ll = 0; ll < n; ll++)
            if (ll != icol) {
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
                for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
            }
    }
    for (l = n - 1; l >= 0; l--) {
        if (indxr[l] != indxc[l])
            for (k = 0; k < n; k++)
                SWAP(a[k][indxr[l]], a[k][indxc[l]]);
    }
}

// From Numerical Recipes Ch 2
void gaussj(MatDoub_IO &a) {
    MatDoub b(a.nrows(), 0);
    gaussj(a, b);
}

void write1(string name, int suffix, function<double(int)> g, int max) {
    stringstream fname;
    fname << name << suffix << ".dat";
    ofstream of;
    of.open(fname.str().c_str());
    for (int i = 0; i < max; i++) {
        of << g(i) << endl;
    }
    of.close();
}

void write2(string name, int suffix, function<double(int)> g, int max) {
    stringstream fname;
    fname << name << suffix << ".dat";
    ofstream of;
    of.open(fname.str().c_str());
    for (int i = 0; i < max; i++) {
        of << ts[i] << " " << g(i) << endl;
    }
    of.close();
}

void write_all_fit_data(vector<math_func> fs, bool seasonal_var) {
    MatDoub_IO F = get_F(fs);
    MatDoub_IO b = get_b(fs);
    gaussj(F, b);
    // now b contains the solution vector (a), and F is (old F)^-1
    MatDoub &F_inv = F;
    MatDoub &a = b;
    int m = seasonal_var ? fs.size() : 3;
    math_func f = [ = , &a](double t) {
        double sum = 0.0;
        for (int k = 0; k < m; k++) {
            sum += a[k][0] * fs[k](t);
        }
        return sum;
    };
    auto get_params = [&a](int i) {
        return a[i][0];
    };
    write1("params", m, get_params, m);
    auto get_fit = [f](int i) {
        return f(ts[i]);
    };
    write2("fit", m, get_fit, N);
    auto get_error = [F_inv](int i) {
        return F_inv[i][i];
    };
    write1("error", m, get_error, m);
    auto get_chi_sq = [f, m](int j) {
        double sum = 0.0;
        for (int i = 0; i < N; i++) {
            double tmp = (ys[i] - f(ts[i])) / sigmas[i];
            sum += tmp * tmp;
        }
        return sum / (N - m);
    };
    write1("chi_sq", m, get_chi_sq, 1);
}

int main() {
    load_arrays();
    transform(ys.begin(), ys.end(), sigmas.begin(), [](double y) {
        return y * sigma_perc;
    });
    auto get_orig_data = [](int i) {
        return ys[i];
    };
    write2("orig", 0, get_orig_data, N);
    // with first harmonic
    vector<math_func> fs = get_fit_funcs(true);
    write_all_fit_data(fs, true);
    // no first harmonic
    vector<math_func> fs2 = get_fit_funcs(false);
    write_all_fit_data(fs2, true);
    // with seasonal variations removed after fit calculation
    write_all_fit_data(fs, false);
}