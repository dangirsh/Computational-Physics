from pylab import *

data = loadtxt('energy.dat', 'float')
x = data[:, 0]
y = data[:, 1]
plot(x, y, 'k-')
xlabel('step')
ylabel('energy')
savefig('energy.png')
