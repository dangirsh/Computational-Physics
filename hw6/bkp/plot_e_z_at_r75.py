from pylab import *

data = loadtxt('e_z_at_r75.dat', 'float')
z = data[:, 0]
e_z = data[:, 1]
plot(z, e_z, 'k-')
xlabel("z (mm)")
ylabel("E_z (V/mm)")
savefig('e_z_at_r75.png')
