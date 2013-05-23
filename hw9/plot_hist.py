from pylab import *

xs = loadtxt('hist.dat', 'float')
hist(xs, 100, normed=True, color = 'w')
xlabel('value ranges')
ylabel('probability')
savefig('hist.png')
