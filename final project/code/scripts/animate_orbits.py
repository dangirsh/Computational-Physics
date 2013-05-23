import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

orbit_name = sys.argv[1]

fig = plt.figure()
ax = fig.add_subplot(111)

data = np.loadtxt('../output/' + orbit_name + '/orbits.dat', 'float')
step_times = data[:, 0]

xs, ys = range(3), range(3)
for i in range(3):
    xs[i] = data[:, i + 1]
    ys[i] = data[:, i + 4]

lines = range(3)
for i in range(3):
    line, = ax.plot(xs[i][:1], ys[i][:1])
    lines[i] = line

time = 0.0  # s
time_index = 0
next_step_time = step_times[0]

interval = 25.0  # ms
last_step_time = step_times[-1]
num_steps = int(last_step_time / (interval / 1000))


def animate(i):
    global time, next_step_time, time_index
    if time > last_step_time:
        time = 0.0  # s
        time_index = 0
        next_step_time = step_times[0]
    while time > next_step_time:
        time_index += 1
        for i in range(3):
            lines[i].set_xdata(xs[i][:time_index])
            lines[i].set_ydata(ys[i][:time_index])
        next_step_time = step_times[time_index]
    time += interval / 1000
    return tuple(lines)


#Init only required for blitting to give a clean slate.
def init():
    for i in range(3):
        line, = ax.plot(xs[i][:1], ys[i][:1])
        lines[i] = line
    ax.set_xbound(-1, 1)
    ax.set_ybound(-1, 1)
    return tuple(lines)

ani = animation.FuncAnimation(fig, animate, np.arange(1, num_steps), init_func=init,
                              interval=interval, blit=True)
# ani.save('../output/' + orbit_name + '/animation.mp4') # TODO: fix
fig.canvas.set_window_title(orbit_name)
plt.show()
