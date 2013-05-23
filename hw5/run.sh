#!/bin/bash

g++ -std=c++0x -O2 -o test hw6.cpp
./test
python countour_plot.py