from pylab import *
import sys
if len(sys.argv) != 5:
    print 'Usage: %s <data file name> <output file name> <body index> <x or y>' % sys.argv[0]
    sys.exit(1)
data_file_name = sys.argv[1]
output_file = sys.argv[2]
body_index = int(sys.argv[3]) + 1
var_name = sys.argv[4]
var_index = 0
if var_name == 'x':
    var_index = 0
elif var_name == 'y':
    var_index = 1
else:
    raise
data = loadtxt(data_file_name, 'float')
plot(data[:, 0], data[:, body_index + 3 * var_index], 'k-')
xlabel("t")
ylabel(var_name)
savefig(output_file)