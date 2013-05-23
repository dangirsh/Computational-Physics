from pylab import *

data = loadtxt('orig.dat', 'float')
t = data[:, 0]
y = data[:, 1]
plot(t, y, 'k-')
xlabel('time (months)')
ylabel('C02 (ppm)')
savefig('orig.png')
