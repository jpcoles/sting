import sys
from pylab import figure, plot, savefig, ion, scatter, draw, xlim, ylim, axis
from time import sleep

ion()

figure()
axis('equal')

X = []
Y = []
l, = plot(X,Y, 'k,')
S = 1
xlim(-1.5*S,1.5*S)
ylim(-1.5*S,1.5*S)

while True:
    p = sys.stdin.readline()
    if not p: break

    p = p.split()
    if not p:
        l.set_data(X,Y)
        draw()
        sleep(0.001)
#        X,Y = [],[]
    else:
        try:
            X.append(float(p[1]))
            Y.append(float(p[2]))
        except ValueError:
            pass

