import sys
from itertools import izip
from pylab import figure, plot, savefig, ion, scatter, draw, xlim, ylim, axis, show, semilogy
from matplotlib.collections import LineCollection
from time import sleep
from numpy import abs, array, sqrt, pi, arange

ion()

    



#axis('equal')
#ax1 = f0.add_subplot(331)
#ax2 = f0.add_subplot(312)
#ax3 = f0.add_subplot(332)
#ax4 = f0.add_subplot(333)
#ax5 = f0.add_subplot(334)


f0 = figure(figsize=(10,10))
if len(sys.argv) > 1:
    f0.suptitle(sys.argv[1])

f1 = figure(figsize=(20,20))
ax1 = f1.add_subplot(111)

ax3 = f0.add_subplot(532)
ax7 = f0.add_subplot(533)

ax2 = f0.add_subplot(512); ax2.set_ylabel('Energy')
ax5 = f0.add_subplot(513); ax5.set_ylabel('Mag. Field')
ax6 = f0.add_subplot(514)
ax8 = f0.add_subplot(515)

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

l1, = ax1.plot(X,Y, 'b,', alpha=1)
l1b, = ax1.plot(X,Y, 'r,', alpha=1)
print l1
S = 10
ax1.axis('scaled')
ax1.set_xlim(-1.5*S,1.5*S)
ax1.set_ylim(-1.5*S,1.5*S)

E = []
J = []
M = []

G = []
K = []
C = []
D = []

Kd = []
Kc = []

# Eccentricities (initial / final)
ei = []
ef = []


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
#ax5.set_ylim(-3.0,3.0)

l4G, = ax6.plot(G, 'k-', lw=2)
l4K, = ax6.plot(K, 'r-', lw=2)
l4C, = ax6.plot(C, 'g-', lw=2)
l4D, = ax6.plot(D, 'b-', lw=2)
ax6.set_xlim(0,3000)
ax6.set_ylim(-2,2)
ax6.set_ylim(-.02,.02)

l5, = ax7.plot([], ef, 'k,')
ax7.set_ylim(-0.5,2)

l6a, = ax8.plot([],Kc, 'r,')
l6b, = ax8.plot([],Kd, 'b,')
ax8.set_xlim(0,3000)
ax8.set_ylim(0,0.11)

T = []
step = 0
frame=0
while True:
    p = sys.stdin.readline()
    if not p: break

    p = p.split()

    if p[0] == 'STEP':
        step = p[1] = int(p[1])
        N    = p[2] = int(p[2])
        p[3:] = map(float, p[3:])

        N    = p[2]
        T.append(p[1])
        E.append(p[3])

        if len(J) and J[0] == 0 and p[8] != 0:
            J[0] = p[8]
        J.append(p[4])
        M.append(p[5])
        G.append(p[6])
        K.append(p[7])
        C.append(p[8])
        D.append(p[9])

        Kc.append(p[10])
        Kd.append(p[11])


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
            l5.set_data(R, ef)
            #l5.set_data(arange(N), ef)
            ax7.set_xlim(0,2)
            #ax7.set_xlim(-1,len(ef))

            l6a.set_data(T,Kc)
            l6b.set_data(T,Kd)

            if J[0] != 0:
                #l2b.set_data(T,map(lambda x: x/J[0], J))
                l2b.set_data(T,map(lambda x: x, J))

            l4G.set_data(T,map(lambda x: x/E[0], G))
            l4K.set_data(T,map(lambda x: x/E[0], K))
            l4C.set_data(T,map(lambda x: x/E[0], C))
            l4D.set_data(T,map(lambda x: x/E[0], D))

            #if M[0] != 0:
                #l3.set_data(T,map(lambda x: x/M[0], M))
            #ax2.set_ylim(E[0]*3, abs(E[0]))
        ax3.clear(); 
        #ax3.hist(R, bins=10, range=(0,2), color='k', histtype='step')
        #if Rm: ax3.hist(Rm, bins=10, range=(0,2), color='b', histtype='step')
        #if Rp: ax3.hist(Rp, bins=10, range=(0,2), color='r', histtype='step')
        #ax3.hist(Xm, bins=10, color='b', histtype='step')
        #ax3.hist(Xp, bins=10, color='r', histtype='step')

        #ax4.clear(); ax4.hist(Y, bins=10, color='k', histtype='step')
        #ax4.hist(Ym, bins=10, color='b', histtype='step')
        #ax4.hist(Yp, bins=10, color='r', histtype='step')
        #print frame
        #savefig('st%05i.png' % frame)
        #frame += 1
        if step % 10 == 0: draw()
        #sleep(0.001)
        #X,Y = [],[]
        #Xm,Ym = [],[]
        #Xp,Yp = [],[]

        tlen = 100
        X,Y = X[-tlen:], Y[-tlen:]
        Xm,Ym = Xm[-tlen:], Ym[-tlen:]
        Xp,Yp = Xp[-tlen:], Yp[-tlen:]
        R = []
        Rm = []
        Rp = []
        ef = []

    for i in xrange(N):
        p = sys.stdin.readline()
        if not p: break

        p = p.split()

        p[0] = int(p[0])
        p[1:] = map(float, p[1:])

        charge = p[5]

        r = sqrt(p[1]**2 + p[2]**2)
        R.append(r)
        X.append(p[1])
        Y.append(p[2])

        ef.append(p[6])

        if charge < 0:
            Xm.append(p[1])
            Ym.append(p[2])
            Rm.append(r)
        else:
            Xp.append(p[1])
            Yp.append(p[2])
            Rp.append(r)


    if not ei:
        ei = ef
