import minimisation
import numpy as np
import os
import analyse
import math
import calculategrad
import initialise
import matplotlib.pyplot as pl


multconst = 1e0
zcomp = 2
Lx,Lz = 35,35
Lt = 100
ws = 0
maxfuncls = 10000
gconv = 2.2e-16
fconv = 2.220446049250313e-16
iters = 1e3
funciters = 1e8
lsiters = 1e9
splitt = 3.0
splitz = 3.0
timesplit = int(Lt / splitt)
timesplitf = int((splitt-1)*(Lt / splitt))
lowz = int(Lz/splitz)
highz = int((splitz-1))*int(Lz/splitz)
dz = 1.0
coza,cozb = 0,0
#time energy
sig = 1e5
#Elastic constants
kt = 1e-15*multconst
ks = 1e-15*multconst
#Bulk constants
a=3e-1
b=2e-1
c=1e-1

#initial and final windings for initial conditions

w0 = 1.0
w1 = 3.0
#Chirality
wideal=(w0+w1)/2.0
#wideal = 1000.0
q0 = wideal*math.pi/(dz*(Lz-1))


#Create data folder
if not os.path.exists('results/ks%ekt%eabc%e%e%esig%dq0%f' % (ks,kt,a,b,c,sig,wideal)):
    os.makedirs('results/ks%ekt%eabc%e%e%esig%dq0%f' % (ks,kt,a,b,c,sig,wideal))

os.chdir('results/ks%ekt%eabc%e%e%esig%dq0%f' % (ks,kt,a,b,c,sig,wideal))
s = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)

Q3 = np.zeros([Lx,Lz,Lt])
Q5 = np.zeros([Lx,Lz,Lt])
Q3temp = np.zeros([Lx,Lz,Lt])
Q5temp = np.zeros([Lx,Lz,Lt])
Q1 = np.zeros([Lx,Lz,Lt])
Q2 = np.zeros([Lx,Lz,Lt])
Q4 = np.zeros([Lx,Lz,Lt])
evalp, evaln = [],[]
i = 0
len = []

minimisation.LBFGS(zcomp, Lx, Lz, Lt, ws, maxfuncls, gconv, fconv, iters, funciters, lsiters, splitt, splitz, dz, coza, cozb, sig, ks, kt, a, b, c, s, w0, w1,multconst,wideal)

#(Q1,Q2,Q3,Q4,Q5) = initialise.initialise(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1,Q2,Q3,Q4,Q5,Lx,Lz,Lt,w0,w1,a,b,c)

data = np.loadtxt('energyarray.dat')

t = 2
for x in range(0,Lx):
    for z in range(0, Lz):
        Q1[x,z,t] = data[5*x + 5*Lz + 5*t*Lx*Lz]
        Q2[x,z,t] = data[5*x + 1 + 5*t*Lx*Lz]
        Q3[x,z,t] = data[5*x + 2 + 5*t*Lx*Lz]
        Q4[x,z,t] = data[5*x + 3 + 5*t*Lx*Lz]
        Q5[x,z,t] = data[5*x + 4 + 5*t*Lx*Lz]

print('Calculating Hessian')
calculategrad.Hessian(sig,Lx,Lz,Lt,ks,kt,q0,s,a,b,c,Q1,Q2,Q3,Q4,Q5)

def teststate():
    for multi in range(-12,12):
        ks = 1*10**(multi)
        if multi % 10 == 2:
            print(multi)
        kt = 1*10**(multi)

        minimisation.LBFGS(zcomp, Lx, Lz, Lt, ws, maxfuncls, gconv, fconv, iters, funciters, lsiters, splitt, splitz, dz, coza, cozb, sig, ks, kt, a, b, c, s, w0, w1,multconst,wideal)
        pos, neg = calculategrad.Hessian(sig,Lx,Lz,Lt,ks,kt,q0,s,a,b,c,Q1,Q2,Q3,Q4,Q5)
        evalp.append(pos)
        evaln.append(neg)
        len.append(multi)
        #c = 1*10**(multi)


    pl.plot(len,evalp, label = 'Positive')
    pl.plot(len,evaln, label = 'Negative')
    pl.legend(loc="upper right")
    pl.xlabel('$k_s, k_t = 1^{x}$')
    pl.ylabel('Number of eigenvalues of Hessian')
    pl.show()
    pl.savefig('hessianevals.pdf')

#minimisation.LBFGS(zcomp, Lx, Lz, Lt, ws, maxfuncls, gconv, fconv, iters, funciters, lsiters, splitt, splitz, dz, coza, cozb, sig, ks, kt, a, b, c, s, w0, w1,multconst,wideal)
