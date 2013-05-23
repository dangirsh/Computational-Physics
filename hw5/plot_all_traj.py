from pylab import *
import sys
if len(sys.argv) != 3:
    print 'Usage: %s <data file name> <output file name>' % sys.argv[0]
    sys.exit(1)
data_file_name = sys.argv[1]
output_file = sys.argv[2]
data = loadtxt(data_file_name, 'float')
p1, = plot(data[:, 1], data[:, 4], 'b-')
p2, = plot(data[:, 2], data[:, 5], 'g-')
p3, = plot(data[:, 3], data[:, 6], 'r-')
xlabel("x(t)")
ylabel("y(t)")
legend([p1, p2, p3], ["Earth", "Moon", "Moon 2"])
savefig(output_file)