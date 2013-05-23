#!/bin/bash

g++ -std=c++0x -O3 -o hw8 hw8.cpp
./hw8

python plot_orig.py
open orig.png

python plot_fits.py
open fit3.png fit5.png fit7.png

python plot_resid.py
open resid.png

python plot_fit_on_data.py
open fit_on_data.png