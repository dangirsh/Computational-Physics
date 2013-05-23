import sys
import numpy as np
import matplotlib.pyplot as plt

orbit_name = sys.argv[1]
fig = plt.figure()
ax = fig.add_subplot(111)
fname = '../output/' + orbit_name + '/fft_x1' + '.dat'
data = np.loadtxt(fname, 'float')
n = len(data)
data = data[:100]
freqs = range(0, 100)
plt.plot(freqs, data, 'k-')
fig.canvas.set_window_title(orbit_name)
plt.show()
fig.savefig(fname.replace('.dat', '.png'))
