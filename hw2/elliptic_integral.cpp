/*
AEP 4380 HW #2
Dan Girshovich
2/7/13
Generates data for elliptic integral with the trapezoid rule and Simpson's rule.
Compile with: g++ -std=c++0x -O2 -o gen_data elliptic_integral.cpp
Tested on OSX 10.8.2 with a Intel Core 2 Duo
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <functional>

using namespace std;

// function defs
double feval(double, double);
typedef std::function<double(double)> integrand_t;
integrand_t get_integrand(double);
double trapezoid_integrator(integrand_t, double, double, int);
double simpson_integrator(integrand_t, double, double, int);
typedef double (*integrator_t)(integrand_t, double, double, int);

void gen_table(integrator_t, double, const char *);
void gen_graph_data(integrator_t, double, double, int, int, const char *);
void output(FILE *, int, double, double);
FILE *open_file_w(const char *);

const double pi = 4.0 * atan( 1.0 );

int main()
{
    // create data for tables comparing the two integrators
    // at different values for x and n
    gen_table(trapezoid_integrator, 0.5, "trap_table1.dat");
    gen_table(trapezoid_integrator, 0.9999, "trap_table2.dat");
    gen_table(simpson_integrator, 0.5, "simp_table1.dat");
    gen_table(simpson_integrator, 0.9999, "simp_table2.dat");

    // Create data for graphs of K(x) using 100 values of x
    // between 0 and 0.9999, and n = 512 for 5 sig figs of accuracy
    int n = 512;
    gen_graph_data(trapezoid_integrator, 0.0, 0.9999, 100, n, "trap_graph.dat");
    gen_graph_data(simpson_integrator, 0.0, 0.9999, 100, n, "simp_graph.dat");
    return (EXIT_SUCCESS);
}

// writes the results of various n values to fname using the integrator on the
// integrand determined by x
void gen_table(integrator_t integrator, double x, const char *fname)
{
    FILE *fp = open_file_w(fname);
    integrand_t integrand = get_integrand(x);
    double integral_val;
    for (int n = 4; n <= 4096; n *= 2)
    {
        integral_val = integrator(integrand, 0, pi / 2, n);
        output(fp, n, x, integral_val);
    }
    fclose(fp);
}

// writes K values to fname by running the integrator over integrands produced
// by num_points values of x between x_min and x_max.
void gen_graph_data(integrator_t integrator, double x_min, double x_max,
                     int num_points, int n, const char *fname)
{
    FILE *fp = open_file_w(fname);
    double dx = (x_max - x_min) / num_points;
    double x, k;
    integrand_t integrand;
    for (int i = 0; i < num_points; i++)
    {
        x = x_min + i * dx;
        integrand = get_integrand(x);
        k = integrator(integrand, 0, pi / 2, n);
        output(fp, n, x, k);
    }
    fclose(fp);
}

// integrates the given integrand from t_min to t_max using n trapezoids
double trapezoid_integrator(integrand_t f, double t_min, double t_max, int n)
{
    double dx = (t_max - t_min) / n;
    double sum = 0.0;
    double t1, t2;
    for (int i = 0; i < n; i++)
    {
        t1 = t_min + i * dx;
        t2 = t_min + (i + 1) * dx;
        sum += 0.5 * dx * (f(t1) + f(t2));
    }
    return sum;
}

// integrates the given integrand from t_min to t_max using n steps of
// Simpson's Rule.
double simpson_integrator(integrand_t f, double t_min, double t_max, int n)
{
    double dx = (t_max - t_min) / n;
    double sum = 0.0;
    double t1, t2, t3;
    for (int i = 0; i < n - 1; i += 2)
    {
        t1 = t_min + i * dx;
        t2 = t_min + (i + 1) * dx;
        t3 = t_min + (i + 2) * dx;
        sum += (1.0 / 3) * dx * (f(t1) + 4 * f(t2) + f(t3));
    }
    return sum;
}

// creates an integrand function for a given value of x
integrand_t get_integrand(double x)
{
    return bind(feval, x, placeholders::_1);
}

// evauluates the function in the elliptic integral
double feval(double x, double t)
{
    return (1.0 / sqrt(1.0 - x * sin(t) * sin(t)));
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

void output(FILE *fp, int n, double x, double k)
{
    // data file for Scilab, Matlab, kaleidaGraph and gnuplot
    fprintf(fp, "%i %16.8f %16.8f \n", n, x, k);
    // printf("%i %f %f \n", n, x, k);
}