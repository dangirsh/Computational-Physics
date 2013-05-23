#!/bin/bash

g++ -std=c++0x -O3 -o hw7 hw7.cpp
./hw7

python plot_psis.py
python plot_potential.py

open psi0_1.png psi1_1.png psi2_1.png psi3_1.png psi4_1.png psi5_1.png
open psi0_2.png psi1_2.png psi2_2.png psi3_2.png psi4_2.png psi5_2.png
open psi0_3.png psi1_3.png psi2_3.png psi3_3.png psi4_3.png psi5_3.png

open potential.png