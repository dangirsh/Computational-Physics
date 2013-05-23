from pylab import *

data = loadtxt('potential.dat', 'float')
x = data[:, 0]
pot = data[:, 1]
plot(x, pot, 'k-')
xlabel('x (A)')
ylabel('potential (V)')
savefig('potential.png')
