from pylab import *
clf()
data = loadtxt("v_contour.dat", 'float')
r = loadtxt("r_contour.dat")
z = loadtxt("z_contour.dat")
cs = contour(z, r, data, 10)  # only lines
clabel(cs)  # label contour levels
xlabel("z (in mm.)")
ylabel("r (in mm.)")
savefig('contour.png')
