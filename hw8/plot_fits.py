from pylab import *

for m in [3, 5, 7]:
    data = loadtxt('fit' + str(m) + '.dat', 'float')
    t = data[:, 0]
    y = data[:, 1]
    clf()
    plot(t, y, 'k-')
    xlabel('time (months)')
    ylabel('C02 (ppm)')
    savefig('fit' + str(m) + '.png')
