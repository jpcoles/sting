import sys
from itertools import izip
from pylab import figure, plot, savefig, ion, scatter, draw, xlim, ylim, axis, show, semilogy
from matplotlib.collections import LineCollection
from time import sleep
from numpy import abs, array, sqrt, pi

ion()

f0 = figure(figsize=(10,10))
#axis('equal')
#ax1 = f0.add_subplot(331)
#ax2 = f0.add_subplot(312)
#ax3 = f0.add_subplot(332)
#ax4 = f0.add_subplot(333)
#ax5 = f0.add_subplot(334)

ax1 = f0.add_subplot(431)
ax2 = f0.add_subplot(412)
ax3 = f0.add_subplot(432)
ax5 = f0.add_subplot(413)
ax6 = f0.add_subplot(414)

X = []
Y = []
C = []
R = []

Xp = []
Yp = []
Xm = []
Ym = []

Rp = []
Rm = []

l1, = ax1.plot(X,Y, 'b,')
l1b, = ax1.plot(X,Y, 'r,')
print l1
S = 1
ax1.set_xlim(-1.5*S,1.5*S)
ax1.set_ylim(-1.5*S,1.5*S)

E = []
J = []
M = []

G = []
K = []
C = []
D = []


l2, = ax2.plot(E, 'k-', lw=2)
l2b, = ax2.plot(J, 'r-', lw=2)
ax2.axhline(1, c='grey', zorder=-1000)
ax2.set_xlim(0, 3000)
#ax2.set_ylim(1e-4, 1e+4)
#ax2.set_ylim(0.90, 1.10)
ax2.set_ylim(0.99, 1.01)

l3, = ax5.plot(M, 'k-', lw=2)
ax5.set_xlim(0,3000)
ax5.set_ylim(-.03,.03)
ax5.set_ylim(-3.0,3.0)

l4G, = ax6.plot(G, 'k-', lw=2)
l4K, = ax6.plot(K, 'r-', lw=2)
l4C, = ax6.plot(C, 'g-', lw=2)
l4D, = ax6.plot(D, 'b-', lw=2)
ax6.set_xlim(0,3000)
ax6.set_ylim(-2,2)
ax6.set_ylim(-.02,.02)

T = []
step = 0
frame=0
while True:
    p = sys.stdin.readline()
    if not p: break

    p = p.split()
    if not p:
        l1.set_data(Xm,Ym)
        l1b.set_data(Xp,Yp)
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
            l3.set_data(T,M)

            if J[0] != 0:
                l2b.set_data(T,map(lambda x: x/J[0], J))

            l4G.set_data(T,map(lambda x: x/E[0], G))
            l4K.set_data(T,map(lambda x: x/E[0], K))
            l4C.set_data(T,map(lambda x: x/E[0], C))
            l4D.set_data(T,map(lambda x: x/E[0], D))

            #if M[0] != 0:
                #l3.set_data(T,map(lambda x: x/M[0], M))
            #ax2.set_ylim(E[0]*3, abs(E[0]))
        ax3.clear(); 
        #ax3.hist(R, bins=10, range=(0,2), color='k', histtype='step')
        ax3.hist(Rm, bins=10, range=(0,2), color='b', histtype='step')
        #ax3.hist(Rp, bins=10, range=(0,2), color='r', histtype='step')
        #ax3.hist(Xm, bins=10, color='b', histtype='step')
        #ax3.hist(Xp, bins=10, color='r', histtype='step')

        #ax4.clear(); ax4.hist(Y, bins=10, color='k', histtype='step')
        #ax4.hist(Ym, bins=10, color='b', histtype='step')
        #ax4.hist(Yp, bins=10, color='r', histtype='step')
        #print frame
        #savefig('st%05i.png' % frame)
        #frame += 1
        draw()
        #sleep(0.001)
        #X,Y = [],[]
        #Xm,Ym = [],[]
        #Xp,Yp = [],[]
    else:
        try:
            p[0] = int(p[0])
            p[1] = int(p[1])
            p[2:] = map(float, p[2:])
            charge = p[6]
            Mfield = p[9]
            if p[0] != step:
                tlen = 1000
                X,Y = X[-tlen:], Y[-tlen:]
                Xm,Ym = Xm[-tlen:], Ym[-tlen:]
                Xp,Yp = Xp[-tlen:], Yp[-tlen:]
                step = int(p[0])
                E.append(p[7])
                if len(J) and J[0] == 0 and p[8] != 0:
                    J[0] = p[8]
                #if len(M) and M[0] == 0 and p[10] != 0:
                    #M[0] = p[8]
                J.append(p[8])
                T.append(step)
                M.append(Mfield)
                R = []
                Rm = []
                Rp = []
                G.append(p[10])
                K.append(p[11])
                C.append(p[12])
                D.append(p[13])

            #r = 2*pi*sqrt(p[2]**2 + p[3]**2)
            r = sqrt(p[2]**2 + p[3]**2)
            R.append(r)
            X.append(p[2])
            Y.append(p[3])


            if charge < 0:
                Xm.append(p[2])
                Ym.append(p[3])
                Rm.append(r)
            else:
                Xp.append(p[2])
                Yp.append(p[3])
                Rp.append(r)

            #C.append(float(p[8])/100)

            #if int(p[1]) == 0:

        except ValueError:
            pass


