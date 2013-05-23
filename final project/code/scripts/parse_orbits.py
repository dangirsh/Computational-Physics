import sys

fname = '../orbits/' + sys.argv[1] + '.orb'
with open(fname, 'r') as f:
    data = f.readlines()

data = filter(lambda l: bool(l) and l != '\n', data)  # remove whitespace

f = lambda i: eval(data[i])

glob_d = f(0)
n = f(1)
period = f(2)

masses = range(n)
positions = range(n)
velocities = range(n)

for i in range(n):
    j = 3 * (i + 1)
    masses[i] = f(j)
    positions[i] = f(j + 1)
    velocities[i] = f(j + 2)

xs, ys = zip(*positions)
vxs, vys = zip(*velocities)

args = [n, period] + masses + list(xs) + list(ys) + list(vxs) + list(vys)

print ' '.join(map(str, args))
