/*
AEP 4380 HW #9
Dan Girshovich
4/18/13
Monte Carlo Calculations
Compile with: g++ -std=c++0x -O2 -o gen_data hw9.cpp
Tested on Mac OSX 10.8.2 with a Intel Core 2 Duo
*/

#include "chain.cpp"

// performs num_steps + 1 monte carlo steps
void monte_carlo(int num_steps) {
    Chain chain;
    ofstream of1, of2;
    of1.open("energy.dat");
    of2.open("length.dat");
    for (int i = 0; i <= num_steps; i++) {
        Chain new_chain(chain);
        new_chain.fold();
        double dE = new_chain.energy - chain.energy;
        if (dE <= 0 || rnd.doub() < exp(-1 * dE)) {
            chain = new_chain;
        }
        if ((i % ((num_steps / 500) + 1)) == 0) {
            of1 << i << " " << chain.energy << std::endl;
            of2 << i << " " << chain.get_length() << std::endl;
            std::cout << "\r" << (double) i / num_steps << std::flush;
        }
        if (i == 1e4 || i == 1e5 || i == 1e6 || i == 1e7) {
            chain.write((int) log10(i));
        }
    }
    of1.close();
    of2.close();
}

void gen_hist_data() {
    ofstream of;
    of.open("hist.dat");
    for (int i = 0; i < 100000; i++) {
        of << rnd.doub() << endl;
    }
    of.close();
}

void gen_points_data() {
    ofstream of;
    of.open("points.dat");
    for (int i = 0; i < 10000; i++) {
        of << rnd.doub() << " " << rnd.doub() << endl;
    }
    of.close();
}

int main() {
    gen_hist_data();
    gen_points_data();
    monte_carlo(1e7);
}