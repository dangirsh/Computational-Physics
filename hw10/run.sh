#!/bin/bash

g++ -I/usr/local/fftw/include/ \
-L/usr/local/fftw/lib/ \
-lfftw3 \
-std=c++0x \
-O2 -Wall -o hw10 hw10.cpp

./hw10

python plot_psis.py
open img/psi_0.png img/psi_0_1.png img/psi_0_2.png img/psi_0_3.png