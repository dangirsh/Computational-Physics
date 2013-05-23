from pylab import *

data = loadtxt('v_at_r75.dat', 'float')
z = data[:, 0]
v = data[:, 1]
plot(z, v, 'k-')
xlabel("z (mm)")
ylabel("voltage (V)")
savefig('v_at_r75.png')
