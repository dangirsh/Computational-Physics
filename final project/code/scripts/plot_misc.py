import numpy as np
import matplotlib.pyplot as plt
import sys

orbit_name = sys.argv[1]
fig = plt.figure()
ax = fig.add_subplot(111)
data = np.loadtxt('../output/' + orbit_name + '/orbits.dat', 'float')
times = data[:, 0]
hs = range(0, len(times) - 2)

for i in range(1, len(times) - 1):
    hs[i - 1] = times[i] - times[i - 1]

plt.plot(times[:len(times) - 2], hs, 'ko')
plt.plot(times[:len(times) - 2], hs, 'k-')
fig.savefig('../output/' + orbit_name + '/h.png')

# plt.clf()

# energies = np.loadtxt('../output/' + orbit_name + '/energy.dat', 'float')
# plt.plot(range(len(energies)), energies, 'k-')
# fig.savefig('../output/' + orbit_name + '/energy.png')
