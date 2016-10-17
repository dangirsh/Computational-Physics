import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

orbit_name = sys.argv[1]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
data = np.loadtxt('../output/' + orbit_name + '/sphere.dat', 'float')
xs = data[:, 0]
ys = data[:, 1]
zs = data[:, 2]
ax.plot(xs, ys, zs, 'k-')
fig.canvas.set_window_title(orbit_name)
plt.show()
# fig.savefig('../output/' + orbit_name + '/sphere.png')
