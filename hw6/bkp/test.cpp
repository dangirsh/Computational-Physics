#include <cstdio>
#include <cstdlib>
#define ARRAYT_BOUNDS_CHECK
#include "array.hpp"

float R_MIN = 0, R_MAX = 15, Z_MAX = 15, Z_MIN = -15;

class solver {
private:
    const float h;
public:
    int I_MAX, J_MAX;
    array<float> v; // values of the potential on the grid
    solver(float _h);
    float get_h();
};

// //assumes valid spacing
//assumes valid spacing
solver::solver(float _h) :
    h(_h), I_MAX((Z_MAX - Z_MIN) / h), J_MAX((R_MAX - R_MIN) / h),
    v(I_MAX + 1, J_MAX + 1) {
    printf("%f %f %d %d %f %f %f %f\n", h, _h, I_MAX, J_MAX, Z_MAX, Z_MIN, R_MAX, R_MIN);
    }

float solver::get_h(){ return h; }

int main()
{
    float h = 0.1;
    solver s(h);
    printf("%f %d %d \n", s.get_h(), s.v.n(), s.v.n());
}