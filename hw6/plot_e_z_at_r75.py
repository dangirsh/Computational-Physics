from pylab import *

data = loadtxt('e_z_at_r75.dat', 'float')
z = data[:, 0]
e_z = map(lambda x: x / 100000, data[:, 1])
plot(z, e_z, 'k-')
xlabel("z (mm)")
ylabel("E_z 10^5 (V/m)")
savefig('e_z_at_r75.png')
