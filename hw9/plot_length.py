from pylab import *

data = loadtxt('length.dat', 'float')
x = data[:, 0]
y = data[:, 1]
plot(x, y, 'k-')
xlabel('step')
ylabel('length')
savefig('length.png')
