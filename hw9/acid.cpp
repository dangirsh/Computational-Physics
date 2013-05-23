#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

Ran rnd(time(NULL));

class Acid {
public:
    int x, y, t;
    Acid(int x, int y, int t): x(x), y(y), t(t) {}
    Acid(): x(-1), y(-1), t(-1) {} // for use as array element
    Acid(Acid &a): x(a.x), y(a.y), t(a.t) {}
    bool is_neighbor(Acid a);
    void perturb();
};

// true if the acid and acid a are vertically or horizontally adjacent
bool Acid::is_neighbor(Acid a) {
    return (abs(a.x - x) + abs(a.y - y)) == 1;
}

// randomly returns positive or negative one
int get_offset() {
    if (rnd.doub() < 0.5) return -1;
    return 1;
}

// moves the acid to a positional diagonally adjacent to it
void Acid::perturb() {
    x += get_offset();
    y += get_offset();
}