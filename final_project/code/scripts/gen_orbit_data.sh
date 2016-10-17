#!/bin/bash

# parse orbit file
a=`python parse_orbits.py $1`
mkdir -p ../output/$1
# echo $1 $a
# # compile and run c++
clang++ -I/usr/local/fftw/include/ -L/usr/local/fftw/lib/ -lfftw3 \
-std=c++0x -stdlib=libc++ -O2 -o ../output/$1/main ../cpp/main.cpp
../output/$1/main $1 $a

python plot_orbits.py $1 && open ../output/$1/orbits.png
python plot_masses.py $1 && open ../output/$1/mass.png
# python plot_sphere.py $1 # && open ../output/$1/sphere.png
# python plot_phases.py $1 && open ../output/$1/phases.png
# python plot_ffts.py $1  && open ../output/$1/fft_x1.png
# python animate_orbits.py $1 &
# python plot_misc.py $1 && open ../output/$1/h.png # && open ../output/$1/energy.png

# clean up
# find .. -name *.dat -delete