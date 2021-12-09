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
def LBFGS(zcomp, Lx, Lz, Lt, ws, maxfuncls, gconv, fconv, iters, funciters, lsiters, splitt, splitz, dz, coza, cozb, sig, ks, kt, a, b, c, s, w0, w1,multconst,wideal):

    timesplit = int(Lt / splitt)
    timesplitf = int((splitt-1)*(Lt / splitt))
    lowz = int(Lz/splitz)
    highz = int((splitz-1))*int(Lz/splitz)
    splay,twist,bend,surface,bulk = 0.0,0.0,0.0,0.0,0.0
    timeen = 0.0
    print(s)
    #Chirality
    q0 = wideal*math.pi/(dz*(Lz-1))
    splay,twist,bend,surface,bulk = 0.0,0.0,0.0,0.0,0.0
    timeen = 0.0


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

    calculategrad.Hessian(sig,Lx,Lz,Lt,ks,kt,q0,s,a,b,c,Q1,Q2,Q3,Q4,Q5)

    #For My stuff
    if zcomp == 1:
        (alpha,beta,gamma,Q4,Q5) = initialise.redefine(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1,Q2,Q3,Q4,Q5,Lx,Lz,Lt,w0,w1,a,b,c)
    if zcomp == 2:
        (Q1,Q2,Q3,Q4,Q5) = initialise.redefinerandom(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1,Q2,Q3,Q4,Q5,Lx,Lz,Lt,w0,w1,a,b,c)

    #For Yucens Data
    # if zcomp == 1:
    #     (Q1,Q2,Q3,Q4,Q5) = initialise.redefine2(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1,Q2,Q3,Q4,Q5,Lz,Lt,w0,w1,a,b,c)
    #print(Q2)
    file1 = open('initialguess.dat', 'w')
    gradfile1 = open('graden.dat', 'w')
    for t in range(0,Lt):
        for x in range(0,Lx):
            for z in range(0,Lz):
                file1.write("%f\n" % Q1[x,z,t])
                file1.write("%f\n" % Q2[x,z,t])
                file1.write("%f\n" % Q3[x,z,t])
                file1.write("%f\n" % Q4[x,z,t])
                file1.write("%f\n" % Q5[x,z,t])
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


    np.savetxt("energyarray.dat", minen.x )
    en.writeenergy(guess,original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk)
    analyse.analysis(Lz,Lt,ks,kt,a,b,c,sig)
    en.writeenergy(guess,original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk)
    #analyse.analysis(Lz,Lt,ks,kt,a,b,c,sig)
    analyse.curl(Lz,Lx,Lt)
    print(Lx,Lz,Lt)

    multconst *= 10
