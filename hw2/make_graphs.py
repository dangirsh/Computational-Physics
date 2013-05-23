# uses numpy, scipy, matplotlib packages
from pylab import *
import sys
# read three-col data file x and two y values
if len(sys.argv) != 4:
    print 'Usage: %s <data file name> <output file name> <graph title>' % sys.argv[0]
    sys.exit(1)

data_file_name = sys.argv[1]
output_file = sys.argv[2]
graph_title = sys.argv[3]
data = loadtxt(data_file_name, 'float')
plot(data[:, 1], data[:, 2], 'k-')
xlabel("x")
ylabel("K(x)")
title(graph_title)
savefig(output_file)
show()
