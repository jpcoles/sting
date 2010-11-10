import sys
from itertools import izip
from pylab import figure, plot, savefig, ion, scatter, draw, xlim, ylim, axis, show, semilogy
from matplotlib.collections import LineCollection
from time import sleep
from numpy import abs, array

ion()

f0 = figure(figsize=(10,10))
#axis('equal')
ax1 = f0.add_subplot(221)
ax2 = f0.add_subplot(212)

X = []
Y = []
C = []
R = []

l1, = ax1.plot(X,Y, 'k,')
print l1
S = 10
ax1.set_xlim(-1.5*S,1.5*S)
ax1.set_ylim(-1.5*S,1.5*S)

E = []
l2, = ax2.plot(E, 'k-')
ax2.axhline(1, c='grey', zorder=-1000)
ax2.set_xlim(0, 1000)
#ax2.set_ylim(1e-4, 1e+4)
ax2.set_ylim(0.99, 1.01)


T = []
step = 0
while True:
    p = sys.stdin.readline()
    if not p: break

    p = p.split()
    if not p:
        l1.set_data(X,Y)
        #ax1.scatter(X,Y,c=C,s=0)
        #l1.set_color(C)
        #l1.set_data(R)
        #l1.set_data = LineCollection(R)
        #for l,c in izip(l1,C):
            #l.set_color(c)

        #l1.set_array(array(C))
        #ax1.add_collection(lc)
        if step > 1:
            l2.set_data(T,map(lambda x: x/E[0], E))
            #ax2.set_ylim(E[0]*3, abs(E[0]))
        draw()
        #sleep(0.001)
#        X,Y = [],[]
    else:
        try:
            if int(p[0]) != step:
                X,Y = [], []
                step = int(p[0])
                E.append(float(p[7]))
                T.append(step)

            R.append([[float(p[2]), float(p[3])]])
            X.append(float(p[2]))
            Y.append(float(p[3]))
            #C.append(float(p[8])/100)

            #if int(p[1]) == 0:

        except ValueError:
            pass

