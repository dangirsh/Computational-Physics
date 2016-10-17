import sys

fname = '../orbits/' + sys.argv[1] + '.dat'
with open(fname, 'r') as f:
    data = f.readlines()

glob_d = eval(data[0])
positions = eval(data[1], glob_d)
velocities = eval(data[2], glob_d)
period = data[3]

out = [{}, 3, period]

for i in range(3):
    out += [1, positions[i], velocities[i], '']

newfname = fname.replace('.dat', '.orb')

with open(newfname, 'w') as f:
    f.write('\n'.join(map(str, out)))
