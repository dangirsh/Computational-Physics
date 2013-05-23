from pylab import *

data1 = loadtxt('orig.dat', 'float')
data2 = loadtxt('fit7.dat', 'float')
t = data1[:, 0]
y1 = data1[:, 1]
y2 = data2[:, 1]
rc('text', usetex=True)
plot(t, y2 - y1, 'ko')
xlabel('time (months)')
ylabel('$\Delta$ C02 (ppm)')
savefig('resid.png')
