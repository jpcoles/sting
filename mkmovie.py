import sys
from pylab import figure, plot, savefig, ion, scatter, draw, xlim, ylim, axis
from time import sleep

ion()

figure()
axis('equal')

X = []
Y = []
l, = plot(X,Y, 'k,')
xlim(-1.5,1.5)
ylim(-1.5,1.5)

while True:
    p = sys.stdin.readline()
    if not p: break

    p = p.split()
    if not p:
        l.set_data(X,Y)
        draw()
        sleep(0.01)
#        X,Y = [],[]
    else:
        try:
            X.append(float(p[1]))
            Y.append(float(p[2]))
        except ValueError:
            pass

