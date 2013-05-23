from pylab import *

data = loadtxt('e_r_at_r75.dat', 'float')
z = data[:, 0]
e_r = data[:, 1]
plot(z, e_r, 'k-')
xlabel("z (mm)")
ylabel("E_r (V/mm)")
savefig('e_r_at_r75.png')
