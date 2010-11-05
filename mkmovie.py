import sys
from pylab import figure, plot, savefig, ion, scatter, draw, xlim, ylim, axis
from time import sleep
from numpy import abs

ion()

f0 = figure()
#axis('equal')
ax1 = f0.add_subplot(211)
ax2 = f0.add_subplot(212)

X = []
Y = []
l1, = ax1.plot(X,Y, 'k,')
S = 10
ax1.set_xlim(-1.5*S,1.5*S)
ax1.set_ylim(-1.5*S,1.5*S)

E = []
l2, = ax2.plot(E, 'k,')
ax2.set_xlim(0, 4000)
#ax2.set_ylim(-E[0]*5, 5*E[0])


T = []
t = 0
while True:
    p = sys.stdin.readline()
    if not p: break

    p = p.split()
    if not p:
        l1.set_data(X,Y)
        l2.set_data(T,E)
        ax2.set_ylim(E[0]*3, abs(E[0]))
        draw()
        sleep(0.001)
#        X,Y = [],[]
    else:
        try:
            X.append(float(p[1]))
            Y.append(float(p[2]))

            if int(p[0]) == 0:
                E.append(float(p[6]))
                t += 1
                T.append(t)

        except ValueError:
            pass

