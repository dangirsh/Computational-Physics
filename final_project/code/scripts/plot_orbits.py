import sys
import numpy as np
import matplotlib.pyplot as plt

orbit_name = sys.argv[1]
fig = plt.figure()
ax = fig.add_subplot(111)
data = np.loadtxt('../output/' + orbit_name + '/orbits.dat', 'float')
plt.plot(data[:, 1], data[:, 4], 'g-')
plt.plot(data[:, 2], data[:, 5], 'b-')
plt.plot(data[:, 3], data[:, 6], 'r-')
# ax.set_xbound(-1, 1)
# ax.set_ybound(-1, 1)
fig.savefig('../output/' + orbit_name + '/orbits.png')
