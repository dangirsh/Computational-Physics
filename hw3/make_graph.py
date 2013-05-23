# read two-col data file
# uses numpy, scipy, matplotlib packages
from pylab import *
import sys
if len(sys.argv) != 4:
    print 'Usage: %s <data file name> <output file name> <y_label>' % sys.argv[0]
    sys.exit(1)

data_file_name = sys.argv[1]
output_file = sys.argv[2]
data = loadtxt(data_file_name, 'float')
plot(data[:, 0], data[:, 1], 'k-')
xlabel("x")
ylabel(sys.argv[3])
savefig(output_file)
#show()
