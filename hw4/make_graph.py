# read two-col data file
# uses numpy, scipy, matplotlib packages
from pylab import *
import sys
if len(sys.argv) != 3:
    print 'Usage: %s <data file name> <output file name>' % sys.argv[0]
    sys.exit(1)

data_file_name = sys.argv[1]
output_file = sys.argv[2]
data = loadtxt(data_file_name, 'float')
plot(data[:, 1], data[:, 2], 'k-')
xlabel("x(t)")
ylabel("y(t)")
savefig(output_file)