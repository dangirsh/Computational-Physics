from pylab import *

for t in [0, 0.1, 0.2, 0.3]:
    clf()
    name = 'psi_' + str(t)
    pixa = loadtxt('dat/' + name + '.dat', 'float')
    pixa = pixa[::-1]
    img = imshow(pixa, aspect='equal', extent=(0, 300, 0, 300))
    img.set_cmap('gray')
    colorbar()
    xlabel('x (in cm)')
    ylabel('y (in cm)')
    title('wave at t = ' + str(t) + ' sec')

    savefig('img/' + name.replace('.', '_') + '.png')
