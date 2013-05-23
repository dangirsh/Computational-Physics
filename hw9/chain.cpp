#include "nr3.h"
#include "ran.h"
#define ARRAYT_BOUNDS_CHECK
#include "array.hpp"
#include "acid.cpp"

const int NUM_ACIDS = 40, NUM_TYPES = 20;
const double E_min = -5, E_max = -2;

class Chain {
public:
    double energy;
    array<Acid> arr;
    array<double> E;
    Chain();
    Chain(Chain &c);
    void fold();
    int get_length();
    void write(int i);
private:
    void init_acids();
    void init_energies();
    double calc_energy();
    bool contains(int x, int y);
    bool try_fold();
};

Chain::Chain() : E(NUM_TYPES, NUM_TYPES), arr(NUM_ACIDS) {
    init_acids();
    init_energies();
    energy = calc_energy();
}

Chain::Chain(Chain &c): E(c.E), arr(c.arr), energy(c.energy) {}

// initializes the chain in a line with random types
void Chain::init_acids() {
    for (int i = 0; i < NUM_ACIDS; i++) {
        int t = (int) (rnd.doub() * NUM_TYPES);
        Acid a(i, 0, t); // horizontal line
        arr(i) = a;
    }
}

// randomly sets the symmetric energy matrix
void Chain::init_energies() {
    for (int i = 0; i < NUM_TYPES; i++) {
        for (int j = i; j < NUM_TYPES; j++) {
            E(i, j) = E_min + (E_max - E_min) * rnd.doub();
            E(j, i) = E(i, j);
        }
    }
}

// returns the enery of the chain using the E matrix
double Chain::calc_energy() {
    double sum = 0.0;
    for (int i = 0; i < NUM_ACIDS - 1; i++) {
        Acid a = arr(i);
        for (int j = i + 1; j < NUM_ACIDS; j++) {
            Acid b = arr(j);
            bool bonded = abs(i - j) == 1;
            if (!bonded) {
                if (a.is_neighbor(b)) {
                    sum += E(a.t, b.t);
                }
            }
        }
    }
    return sum;
}

// returns true if the chain contains an acid at point (x, y)
bool Chain::contains(int x, int y) {
    for (int i = 0; i < NUM_ACIDS; i++) {
        Acid a = arr(i);
        if (a.x == x && a.y == y) return true;
    }
    return false;
}

// tries to make a random fold and updates the chain and energy if successful
// returns true if successful, false otherwise
bool Chain::try_fold() {
    // randomly pick an acid
    int i = (int) (NUM_ACIDS * rnd.doub());
    Acid new_a(arr(i));
    // randomly choose a perturbation to try
    new_a.perturb();
    // check availability
    if (contains(new_a.x, new_a.y)) return false;
    // check left neighbor correctly bonded
    if (i > 0) {
        Acid left_n = arr(i - 1);
        if (!new_a.is_neighbor(left_n)) return false;
    }
    // check right neighbor correctly bonded
    if (i < NUM_ACIDS - 1) {
        Acid right_n = arr(i + 1);
        if (!new_a.is_neighbor(right_n)) return false;
    }
    // all checks passed
    arr(i) = new_a;
    energy = calc_energy();
    return true;
}

void Chain::fold() {
    while (!try_fold());
}

int Chain::get_length() {
    int min = 1000, max = -1000;
    for (int i = 0; i < NUM_ACIDS; i++) {
        int x = arr(i).x;
        if (x < min) min = x;
        if (x > max) max = x;
    }
    return max - min;
}

void Chain::write(int i) {
    stringstream fname;
    fname << "chain" << i << ".dat";
    ofstream of;
    of.open(fname.str().c_str());
    for (int i = 0; i < NUM_ACIDS; i++) {
        of << arr(i).x << " " << arr(i).y << std::endl;
    }
    of.close();
}