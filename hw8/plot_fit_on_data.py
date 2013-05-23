from pylab import *

data = loadtxt('fit7.dat', 'float')
data2 = loadtxt('orig.dat', 'float')
t = data[:, 0]
y = data[:, 1]
y2 = data2[:, 1]
clf()
plot(t, y, 'r-')
plot(t, y2, 'k.')
xlabel('time (months)')
ylabel('C02 (ppm)')
savefig('fit_on_data.png')
