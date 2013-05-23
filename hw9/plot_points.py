from pylab import *

data = loadtxt('points.dat', 'float')
x = data[:, 0]
y = data[:, 1]
plot(x, y, 'k,')
xlabel('x')
ylabel('y')
savefig('points.png')
