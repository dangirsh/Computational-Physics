from pylab import *
import sys
if len(sys.argv) != 3:
    print 'Usage: %s <data file name> <output file name>' % sys.argv[0]
    sys.exit(1)
data_file_name = sys.argv[1]
output_file = sys.argv[2]
data = loadtxt(data_file_name, 'float')
times = data[:, 0]
hs = range(0, len(times))
interesting_times = []
for i in range(0, len(times)-1):
    hs[i] = times[i+1] - times[i]
    if hs[i] < 500:
        interesting_times.append(times[i])
hs[len(times) - 1] = hs[len(times) - 2]
plot(times, hs , 'k-')
xlabel("t")
ylabel("h")
ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
savefig(output_file)

output_file = "img/interesting_points.png"
xe = []
ye = []
xm = []
ym = []
xm2 = []
ym2 = []
for line in data:
    if line[0] in interesting_times:
        xe.append(line[1])
        ye.append(line[2])
        xm.append(line[3])
        ym.append(line[4])
        xm2.append(line[5])
        ym2.append(line[6])
p1, = plot(xe, ye, 'bo')
p2, = plot(xm, ym, 'go')
p3, = plot(xm2, ym2, 'ro')
xlabel("x(t)")
ylabel("y(t)")
legend([p1, p2, p3], ["Earth", "Moon", "Moon 2"])
savefig(output_file)