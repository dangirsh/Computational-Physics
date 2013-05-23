from pylab import *

data_file_name = 'omega.dat'
data = loadtxt(data_file_name, 'float')
omegas = data[:, 0]
ns = data[:, 1]
plot(omegas, ns, 'k-')
xlabel("omega")
ylabel("# iterations")
savefig('n_omega.png')
