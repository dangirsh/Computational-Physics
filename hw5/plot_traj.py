from pylab import *
from itertools import islice
import sys

if len(sys.argv) != 4:
    print 'Usage: %s <data file name> <output file name> <body index>' % sys.argv[0]
    sys.exit(1)
data_file_name = sys.argv[1]
output_file = sys.argv[2]
body_index = int(sys.argv[3]) + 1
data = loadtxt(data_file_name, 'float')
xs = data[:, body_index]
ys = data[:, body_index + 3]
plot(xs, ys, 'k-')
xlabel("x(t)")
ylabel("y(t)")
savefig(output_file)