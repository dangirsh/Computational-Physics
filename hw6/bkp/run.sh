#!/bin/bash

g++ -std=c++0x -O2 -o test hw6.cpp
./test

#part one
# python plot_n_omega.py
# open n_omega.png

#part three
python plot_v_on_axis.py
python plot_e_z_on_axis.py
open v_on_axis.png e_z_on_axis.png

#part four
python plot_v_at_r75.py
python plot_e_z_at_r75.py
python plot_e_r_at_r75.py
open v_at_r75.png e_z_at_r75.png e_r_at_r75.png

# #part 5
python contour_plot.py
open contour.png