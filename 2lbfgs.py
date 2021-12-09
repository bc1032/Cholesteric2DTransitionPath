#Libraries
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize

import periodic
import calculategrad
import initialise
import update
import energy as en
import analyse
import writeenergy

multconst = 1e0
#ks = 1e-12*multconst
#kt = 1e-12*multconst

#while kt <= 1e5:
zcomp = 0
Lx,Lz=20,20
Lt=80
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
ks = 1e-2*multconst
kt = 1e-2*multconst
#Bulk constants
a=3e4
b=2e4
c=1e4
#Create data folder
if not os.path.exists('results/ks%ekt%eabc%e%e%esig%d' % (ks,kt,a,b,c,sig)):
    os.makedirs('results/ks%ekt%eabc%e%e%esig%d' % (ks,kt,a,b,c,sig))


splay,twist,bend,surface,bulk = 0.0,0.0,0.0,0.0,0.0
timeen = 0.0

s = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
print(s)

#initial and final windings for initial conditions

w0 = 1.0
w1 = 3.0
#Chirality
wideal=(w0+w1)/2.0
q0 = wideal*math.pi/(dz*(Lz-1))

fq5 = np.zeros([Lx,Lz,Lt])
fq3 = np.zeros([Lx,Lz,Lt])
fgamma = np.zeros([Lx,Lz,Lt])
falpha = np.zeros([Lx,Lz,Lt])
fbeta = np.zeros([Lx,Lz,Lt])
alpha = np.zeros([Lx,Lz,Lt])
beta = np.zeros([Lx,Lz,Lt])
gamma = np.zeros([Lx,Lz,Lt])
alphatemp = np.zeros([Lx,Lz,Lt])
betatemp = np.zeros([Lx,Lz,Lt])
gammatemp = np.zeros([Lx,Lz,Lt])
Q3 = np.zeros([Lx,Lz,Lt])
Q5 = np.zeros([Lx,Lz,Lt])
Q3temp = np.zeros([Lx,Lz,Lt])
Q5temp = np.zeros([Lx,Lz,Lt])
Q1 = np.zeros([Lx,Lz,Lt])
Q2 = np.zeros([Lx,Lz,Lt])
Q4 = np.zeros([Lx,Lz,Lt])

(Q1,Q2,Q3,Q4,Q5) = initialise.initialise(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1,Q2,Q3,Q4,Q5,Lx,Lz,Lt,w0,w1,a,b,c)

#For My stuff
if zcomp == 1:
    (alpha,beta,gamma,Q4,Q5) = initialise.redefine(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1,Q2,Q3,Q4,Q5,Lx,Lz,Lt,w0,w1,a,b,c)


#For Yucens Data
# if zcomp == 1:
#     (Q1,Q2,Q3,Q4,Q5) = initialise.redefine2(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1,Q2,Q3,Q4,Q5,Lz,Lt,w0,w1,a,b,c)
#print(Q2)
file1 = open('initialguess.dat', 'w')
gradfile1 = open('graden.dat', 'w')
for t in range(0,Lt):
    for x in range(0,Lx):
        for i in range(0,Lz):
            file1.write("%f\n" % Q1[x,i,t])
            file1.write("%f\n" % Q2[x,i,t])
            file1.write("%f\n" % Q3[x,i,t])
            file1.write("%f\n" % Q4[x,i,t])
            file1.write("%f\n" % Q5[x,i,t])
            gradfile1.write("0.0\n")
            gradfile1.write("0.0\n")
            gradfile1.write("0.0\n")
            gradfile1.write("0.0\n")
            gradfile1.write("0.0\n")

file1.close()
gradfile1.close()
GradE = np.loadtxt("graden.dat")
guess = np.loadtxt("initialguess.dat")
original = guess

print(np.shape(guess), type(guess))
z = 0
x = 0
t = 0
#GradE = calculategrad.calcgrad(guess,original,GradE,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws)

#Q3 = np.zeros([Lx,Lz,Lt])
#Q5 = np.zeros([Lx,Lz,Lt])
#Q1 = np.zeros([Lx,Lz,Lt])
#Q2 = np.zeros([Lx,Lz,Lt])
#Q4 = np.zeros([Lx,Lz,Lt])


minen = scipy.optimize.minimize(en.calcenergy,guess,\
args=(original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk),\
options={'disp': True, 'maxiter': iters},method='L-BFGS-B',jac=True)


np.savetxt("results/ks%ekt%eabc%e%e%esig%d/energyarray.dat" % (ks,kt,a,b,c,sig), minen.x )
en.writeenergy(guess,original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk)
#analyse.analysis(Lz,Lt,ks,kt,a,b,c,sig)
print(Lx,Lz,Lt)


# Energy = np.zeros((Lt))
# s = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
# qt = np.array([[2.0*s/3.0,0,0],[0,-s/3.0,0],[0,0,-s/3.0]])
# Qt1,Qt2,Qt3,Qt4,Qt5= (2.0*s/3.0),0.0,0.0,-s/3.0,0.0
#np.savetxt("differencearray.dat", guess-minen.x)

multconst *= 10
