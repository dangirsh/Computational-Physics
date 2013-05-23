from pylab import *

data = loadtxt('v_on_axis.dat', 'float')
z = data[:, 0]
v = data[:, 1]
plot(z, v, 'k-')
xlabel("z (mm)")
ylabel("voltage (V)")
savefig('v_on_axis.png')
