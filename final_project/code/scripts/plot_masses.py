import sys
import numpy as np
import matplotlib.pyplot as plt

orbit_name = sys.argv[1]
fig = plt.figure()
ax = fig.add_subplot(111)
data = np.loadtxt('../output/' + orbit_name + '/mass.dat', 'float')
x = range(0, len(data))
plt.plot(x, data[:, 0], 'g-')
plt.plot(x, data[:, 1], 'b-')
plt.plot(x, data[:, 2], 'r-')
fig.savefig('../output/' + orbit_name + '/mass.png')
