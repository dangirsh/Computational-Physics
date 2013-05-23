#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "nr3.h"
#include "ran.h"
#define ARRAYT_BOUNDS_CHECK
#include "array.hpp"

class A {
public:
    int x, y, z;
    A(): x(-1), y(-1), z(-1) {}
};

int main() {
    array<A> arr(4);
    A a = arr(0);
    A b = A();
    b.x += 1;
    b.y -= 1;
    b.z = 15;
    arr(0) = b;
    a = b;
    std::cout << a.x << " " << a.y << " " << a.z << std::endl;
}