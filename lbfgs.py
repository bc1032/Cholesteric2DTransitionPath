import minimisation
import numpy as np
import os
import analyse
import math
import calculategrad
import initialise

multconst = 1e0
zcomp = 2
Lx,Lz = 20,20
Lt = 10
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
kt = 1e-2*multconst
ks = 1e-2*multconst
#Bulk constants
a=3e-2
b=2e-2
c=1e-2

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

(Q1,Q2,Q3,Q4,Q5) = initialise.initialise(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1,Q2,Q3,Q4,Q5,Lx,Lz,Lt,w0,w1,a,b,c)

calculategrad.Hessian(sig,Lx,Lz,Lt,ks,kt,q0,s,a,b,c,Q1,Q2,Q3,Q4,Q5)

#minimisation.LBFGS(zcomp, Lx, Lz, Lt, ws, maxfuncls, gconv, fconv, iters, funciters, lsiters, splitt, splitz, dz, coza, cozb, sig, ks, kt, a, b, c, s, w0, w1,multconst,wideal)
#np.savetxt("energyarray.dat", minen.x )
# en.writeenergy(guess,original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk)
# analyse.analysis(Lz,Lt,ks,kt,a,b,c,sig)
# analyse.curl(Lz,Lx,Lt)
print(Lx,Lz,Lt)
