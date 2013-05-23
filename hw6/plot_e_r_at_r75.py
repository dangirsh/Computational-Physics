from pylab import *

data = loadtxt('e_r_at_r75.dat', 'float')
z = data[:, 0]
e_r = map(lambda x: x / 10000, data[:, 1])
plot(z, e_r, 'k-')
xlabel("z (mm)")
ylabel("E_r 10^4 (V/m)")
savefig('e_r_at_r75.png')
