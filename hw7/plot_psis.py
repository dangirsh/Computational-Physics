from pylab import *

names = ["$Re[\Psi(x)]$", "$Im[\Psi(x)]$", "$|\Psi(x)|^2$"]

for j in range(1, 4):
    for i in range(0, 6):
        data = loadtxt('psi' + str(i) + '.dat', 'float')
        x = data[:, 0]
        phi = data[:, j]
        clf()
        rc('text', usetex=True)
        plot(x, phi, 'k-')
        xlabel('x (A)')
        ylabel(names[j - 1])
        savefig('psi' + str(i) + "_" + str(j) + '.png')