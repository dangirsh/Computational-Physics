from pylab import *

for i in [4, 5, 6, 7]:
    name = 'chain' + str(i)
    data = loadtxt(name + '.dat', 'int')
    x = data[:, 0]
    y = data[:, 1]
    clf()
    plot(x, y, 'ko')
    plot(x, y, 'k-')
    xlabel('x')
    ylabel('y')
    savefig(name + '.png')
