import sys
import numpy as np
import matplotlib.pyplot as plt

orbit_name = sys.argv[1]
fig = plt.figure()
ax = fig.add_subplot(111)
data = np.loadtxt('../output/' + orbit_name + '/phases.dat', 'float')
colors = ['g', 'b', 'r', 'c', 'm', 'y']
for i in range(6):
    plt.plot(data[:, i], data[:, i + 6], colors[i] + 'o')
fig.savefig('../output/' + orbit_name + '/phases.png')
