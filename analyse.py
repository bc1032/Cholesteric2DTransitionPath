import numpy as np
import pandas as pd
import pylab as pl
import math
import sys
import os
import energy as en
import matplotlib.pyplot as plt
import math as m
import cmath
from mpl_toolkits.mplot3d import Axes3D

os.chdir('/home/bc1032/Desktop/Work/Cholesterics/2DLBFGS/BlueCrystal')

def analysis(lx,lz,lt,ks,kt,a,b,c,sig):
    frames = 1

    data = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyarray.dat" % (ks,kt,a,b,c,sig))

    if not os.path.exists('ks%ekt%eabc%e%e%esig%d/data' % (ks,kt,a,b,c,sig) ):
        os.makedirs('ks%ekt%eabc%e%e%esig%d/data' % (ks,kt,a,b,c,sig))

    r,theta,phi=[],[],[]
    def cart2sph(x,y,z):
        XsqPlusYsq = x**2 + y**2
        r2 = m.sqrt(XsqPlusYsq + z**2)               # r
        elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
        az = m.atan2(y,x)                          # phi
        r.append(r2)
        theta.append(elev)
        phi.append(az)
        phifile.write("%f\n" % (az))
        rfile.write("%f\n" % r2 )
        thetafile.write("%f\n" % elev )
        sphfile.write("%f   %f  %f\n" % (r2,elev,az))
        #print(r2,elev,az)
        return (r2, elev, az)

    cent = int(lt*lz/2.0)
    def complexe_modulo(z):
         a = z.real
         b = z.imag
         return math.sqrt(a**2+b**2)
    def JMatrix(x,y,z):
        beta = m.atan2(y,z)                       # azimuthal angle
        #XsqPlusYsq = x**2 + y**2
        #beta = m.atan2(z,m.sqrt(XsqPlusYsq))
        #r2 = m.sqrt(XsqPlusYsq + z**2)               # r
        Ezi = 10
        Exi = 10
        #alpha = m.atan2(Exi,Eyi)
        gamma = m.atan2(-x,y) #  + alpha
        lamda = 450-9
        n0 = 1.52
        ne = 1.9
        neb = math.sqrt((ne*n0)/(math.sqrt((n0**2 * math.sin(beta)**2 + ne**2*math.cos(beta)**2))))
        #print(beta,gamma)
        phi0 = (2.0*math.pi/lamda)*n0#/(lz)
        phie = (2.0*math.pi/lamda)*neb#/(lz)
        J1 = math.cos(gamma)**2 * cmath.exp(complex(0,1)*phie) + math.sin(gamma)**2 * cmath.exp(complex(0,1)*phi0)
        J2 = math.cos(gamma)*math.sin(gamma)*(cmath.exp(complex(0,1)*phie) - cmath.exp(complex(0,1)*phi0))
        J3 = J2
        J4 = math.sin(gamma)**2 * cmath.exp(complex(0,1)*phie) + math.cos(gamma)**2 * cmath.exp(complex(0,1)*phi0)
        #print(J4)
        J = np.array([[J1,J2],[J3,J4]])
        #print(J)
        matl = np.array([[0,0],[0,1]])
        matr = np.array([[1,0],[0,0]])
        Ein = np.array([[Exi],[Ezi]])
        Eout = np.linalg.multi_dot([matl,J,matr,Ein])
        #print(Eout)
        I = complexe_modulo(Eout[0,0]) + complexe_modulo(Eout[1,0])
        intensityfile.write("%f\n" % I)
        jmats.write("%f %f  %f  %f\n" % (J[0,0],J[0,1],J[1,0],J[1,1]))
        return(I)

    #data = np.loadtxt("string")
    print("File Loaded")
    mod = []
    nx,ny,nz=[],[],[]

    phifile = open("ks%ekt%eabc%e%e%esig%d/data/phi.dat" % (ks,kt,a,b,c,sig), 'w' )
    rfile = open("ks%ekt%eabc%e%e%esig%d/data/r.dat" % (ks,kt,a,b,c,sig), 'w' )
    thetafile = open("ks%ekt%eabc%e%e%esig%d/data/theta.dat" % (ks,kt,a,b,c,sig), 'w' )
    sphfile = open("ks%ekt%eabc%e%e%esig%d/data/sphericalcoords.dat" % (ks,kt,a,b,c,sig), 'w' )
    cartfile = open("ks%ekt%eabc%e%e%esig%d/data/cartesiancoords.dat" % (ks,kt,a,b,c,sig), 'w' )
    evalfile = open("ks%ekt%eabc%e%e%esig%d/data/eval.dat" % (ks,kt,a,b,c,sig), 'w' )
    threedfile1 = open("ks%ekt%eabc%e%e%esig%d/data/3dfile1.dat" % (ks,kt,a,b,c,sig), 'w' )
    threedfile2 = open("ks%ekt%eabc%e%e%esig%d/data/3dfile2.dat" % (ks,kt,a,b,c,sig), 'w' )
    threedfile3 = open("ks%ekt%eabc%e%e%esig%d/data/3dfile3.dat" % (ks,kt,a,b,c,sig), 'w' )
    sumqfile = open("ks%ekt%eabc%e%e%esig%d/data/sumQ.dat" % (ks,kt,a,b,c,sig), 'w' )
    eval1file = open("ks%ekt%eabc%e%e%esig%d/data/eval1.dat" % (ks,kt,a,b,c,sig), 'w' )
    eval2file = open("ks%ekt%eabc%e%e%esig%d/data/eval2.dat" % (ks,kt,a,b,c,sig), 'w' )
    eval3file = open("ks%ekt%eabc%e%e%esig%d/data/eval3.dat" % (ks,kt,a,b,c,sig), 'w' )
    nxfile = open("ks%ekt%eabc%e%e%esig%d/data/nx.dat" % (ks,kt,a,b,c,sig), 'w' )
    nyfile = open("ks%ekt%eabc%e%e%esig%d/data/ny.dat" % (ks,kt,a,b,c,sig), 'w' )
    nzfile = open("ks%ekt%eabc%e%e%esig%d/data/nz.dat" % (ks,kt,a,b,c,sig), 'w' )
    Q1file = open("ks%ekt%eabc%e%e%esig%d/data/Q1.dat" % (ks,kt,a,b,c,sig), 'w' )
    Q2file = open("ks%ekt%eabc%e%e%esig%d/data/Q2.dat" % (ks,kt,a,b,c,sig), 'w' )
    Q3file = open("ks%ekt%eabc%e%e%esig%d/data/Q3.dat" % (ks,kt,a,b,c,sig), 'w' )
    Q4file = open("ks%ekt%eabc%e%e%esig%d/data/Q4.dat" % (ks,kt,a,b,c,sig), 'w' )
    Q5file = open("ks%ekt%eabc%e%e%esig%d/data/Q5.dat" % (ks,kt,a,b,c,sig), 'w' )
    modqfile = open("ks%ekt%eabc%e%e%esig%d/data/modQ.dat" % (ks,kt,a,b,c,sig), 'w')
    intensityfile = open("ks%ekt%eabc%e%e%esig%d/data/Intensity.dat" % (ks,kt,a,b,c,sig), 'w')
    jmats = open("ks%ekt%eabc%e%e%esig%d/data/Jmats.dat" % (ks,kt,a,b,c,sig), 'w')

    Q1m,Q2m,Q3m,Q4m,Q5m = [],[],[],[],[]
    eval1,eval2,eval3 = [],[],[]

    print(np.shape(data))

    for t in range(0,lt):
        for x in range(0,lx):
            for z in range(0,(lz)):

                Q1 = data[int(5*x*lz + 5*z + 5*t*lx*lz)]
                Q2 = data[int(5*x*lz + 5*z + 5*t*lx*lz + 1)]
                Q3 = data[int(5*x*lz + 5*z + 5*t*lx*lz + 2)]
                Q4 = data[int(5*x*lz + 5*z + 5*t*lx*lz + 3)]
                Q5 = data[int(5*x*lz + 5*z + 5*t*lx*lz + 4)]
                M = np.array([[Q1,Q2,Q3],[Q2,Q4,Q5],[Q3,Q5,-Q1-Q4]])

                mod.append((Q1**2 + Q2**2 + Q3**2 + Q4**2 + Q5**2 + Q1*Q4))
                modqfile.write("%f\n" % (Q1**2 + Q2**2 + Q3**2 + Q4**2 + Q5**2 + Q1*Q4))
                sumqfile.write("%f\n" % (Q1+Q2+Q3+Q4+Q5))

                eval, evec = np.linalg.eig(M)
                idx = eval.argsort()[::-1]
                eval = eval[idx]
                evec = evec[:,idx]
                x = abs(evec[0,0])
                y = abs(evec[1,0])
                z = abs(evec[2,0])
                #print(x,y,z)
                xi = evec[0,0]
                yi = (evec[1,0])
                zi = (evec[2,0])
                #JMatrix(x,y,z)

                nx.append(x)
                ny.append(y)
                nz.append(z)
                eval1file.write("%f\n" % eval[0])
                eval2file.write("%f\n" % eval[1])
                eval3file.write("%f\n" % eval[2])
                Q1m.append(Q1)
                Q2m.append(Q2)
                Q3m.append(Q3)
                Q4m.append(Q4)
                Q5m.append(Q5)
                Q1file.write("%f\n" % Q1)
                Q2file.write("%f\n" % Q2)
                Q3file.write("%f\n" % Q3)
                Q4file.write("%f\n" % Q4)
                Q5file.write("%f\n" % Q5)
                eval1.append(eval[0])
                eval2.append(eval[1])
                eval3.append(eval[2])
                evalfile.write("%f  %f  %f\n" % (eval[0],eval[1],eval[2]))
                nxfile.write("%f\n" % x)
                nyfile.write("%f\n" % y)
                nzfile.write("%f\n" % z)
                cartfile.write("%f  %f  %f\n" % (x,y,z))
                cart2sph(x,y,z)
                threedfile1.write("%f\n" % (Q1-Q4) )
                threedfile2.write("%f\n" % (Q2 ))
                threedfile3.write("%f\n" % (Q1+Q4) )

    thetafile.close()
    phifile.close()
    rfile.close()
    sphfile.close()
    cartfile.close()
    Q1file.close()
    Q2file.close()
    Q3file.close()
    Q4file.close()
    Q5file.close()
    modqfile.close()
    sumqfile.close()
    evalfile.close()
    eval1file.close()
    eval2file.close()
    eval3file.close()
    threedfile1.close()
    threedfile2.close()
    threedfile3.close()
    nxfile.close()
    nyfile.close()
    nzfile.close()
    intensityfile.close()
    jmats.close()

    return(0)


def writeenergy(guess,sig,Lz,Lt,ks,kt,a,b,c):
    #print("Energy Called")
    original = guess
    dz = 1
    dt = 1
    w0 = 1.0
    w1 = 3.0
    ws = 0.0
    s = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
    #Chirality
    wideal=(w0+w1)/2.0
    q0 = wideal*math.pi/(dz*(Lz-1))
    Q3 = np.zeros([Lz,Lt])
    Q5 = np.zeros([Lz,Lt])
    Q1 = np.zeros([Lz,Lt])
    Q2 = np.zeros([Lz,Lt])
    Q4 = np.zeros([Lz,Lt])
    bulk, twist, splay = 0.0,0.0,0.0
    Qt1,Qt2,Qt3,Qt4,Qt5= (2.0*s/3.0),0.0,0.0,-s/3.0,0.0

    enarrayfinal = []

    for t in range(0,Lt):
        for z in range(0,Lz):
            Q1[z,t] = guess[5*t*Lz + 5*z]
            Q2[z,t] = guess[5*t*Lz + 5*z + 1]
            Q3[z,t] = guess[5*t*Lz + 5*z + 2]
            Q4[z,t] = guess[5*t*Lz + 5*z + 3]
            Q5[z,t] = guess[5*t*Lz + 5*z + 4]
    enarrayinit = []
    bulkenarr, splayenarr, twistenarr, timeenarr = [],[],[],[]
    for t in range(2,Lt-2):
        splayen,twisten,benden,surfaceen,bulken, timeen, surface = 0.0,0.0,0.0,0.0,0.0,0.0,0.0
        for z in range(2,Lz-2):
            if t == 0:
                Q1[z,t] = original[5*z]
                Q2[z,t] = original[5*z + 1]
                Q3[z,t] = original[5*z + 2]
                Q4[z,t] = original[5*z + 3]
                Q5[z,t] = original[5*z + 4]
            if t == Lt-1:
                Q1[z,t] = original[5*t*Lz + 5*z]
                Q2[z,t] = original[5*t*Lz + 5*z + 1]
                Q3[z,t] = original[5*t*Lz + 5*z + 2]
                Q4[z,t] = original[5*t*Lz + 5*z + 3]
                Q5[z,t] = original[5*t*Lz + 5*z + 4]

            if t != 0 and t != Lt-1:
                timeen += (sig*((-Q1[z,t-1] + Q1[z,t+1])**2/(4.*dz**2) + (-Q2[z,t-1] + Q2[z,t+1])**2/(2.*dz**2)\
                        + (-Q3[z,t-1] + Q3[z,t+1])**2/(2.*dz**2) + \
                        (-Q4[z,t-1] + Q4[z,t+1])**2/(4.*dz**2) + \
                        (-(-Q1[z,t-1] + Q1[z,t+1])/(2.*dz) - (-Q4[z,t-1] + Q4[z,t+1])/(2.*dz))**2 +\
                        (-Q5[z,t-1] + Q5[z,t+1])**2/(2.*dz**2)))/2.

            if z == 0:
                Q1[z,t] = original[5*t*Lz + 5*z]
                Q2[z,t] = original[5*t*Lz + 5*z + 1]
                Q3[z,t] = original[5*t*Lz + 5*z + 2]
                Q4[z,t] = original[5*t*Lz + 5*z + 3]
                Q5[z,t] = original[5*t*Lz + 5*z + 4]


                surface += 0.5*ws*((Q1[z,t]-Qt1)**2 + 2*(Q2[z,t]-Qt2)**2 \
                + 2*(Q3[z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[z,t] - Q4[z,t])**2 \
                + (Q4[z,t] - Qt4)**2 +2*(Q5[z,t] - Qt5)**2 )

                bulk += -(a*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)) +\
                        c*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)**2 +\
                        b*(Q1[z,t]**2*Q4[z,t] + (-Q2[z,t]**2 + Q3[z,t]**2)*Q4[z,t] - 2*Q2[z,t]*Q3[z,t]*Q5[z,t] +\
                        Q1[z,t]*(-Q2[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2))
            if z == Lz-1:
                Q1[z,t] = original[5*t*Lz + 5*z]
                Q2[z,t] = original[5*t*Lz + 5*z + 1]
                Q3[z,t] = original[5*t*Lz + 5*z + 2]
                Q4[z,t] = original[5*t*Lz + 5*z + 3]
                Q5[z,t] = original[5*t*Lz + 5*z + 4]


                surface += 0.5*ws*( (Q1[z,t]-Qt1)**2 + 2*(Q2[z,t]-Qt2)**2 \
                + 2*(Q3[z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[z,t] - Q4[z,t])**2 \
                + (Q4[z,t] - Qt4)**2 +2*(Q5[z,t] - Qt5)**2 )

                bulk += -(a*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)) +\
                        c*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)**2 +\
                        b*(Q1[z,t]**2*Q4[z,t] + (-Q2[z,t]**2 + Q3[z,t]**2)*Q4[z,t] - 2*Q2[z,t]*Q3[z,t]*Q5[z,t] +\
                        Q1[z,t]*(-Q2[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2))

            else:


                bulk += -(a*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)) +\
                        c*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)**2 +\
                        b*(Q1[z,t]**2*Q4[z,t] + (-Q2[z,t]**2 + Q3[z,t]**2)*Q4[z,t] - 2*Q2[z,t]*Q3[z,t]*Q5[z,t] +\
                        Q1[z,t]*(-Q2[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2))

                bulken += -(a*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)) +\
                        c*(Q1[z,t]**2 + Q2[z,t]**2 + Q3[z,t]**2 + Q1[z,t]*Q4[z,t] + Q4[z,t]**2 + Q5[z,t]**2)**2 +\
                        b*(Q1[z,t]**2*Q4[z,t] + (-Q2[z,t]**2 + Q3[z,t]**2)*Q4[z,t] - 2*Q2[z,t]*Q3[z,t]*Q5[z,t] +\
                        Q1[z,t]*(-Q2[z,t]**2 + Q4[z,t]**2 + Q5[z,t]**2))


                twist += (kt*((-Q1[z-1,t] + Q1[z+1,t])**2/(4.*dt**2) - \
                    (2*q0*Q1[z,t]*(-Q2[z-1,t] + Q2[z+1,t]))/dt + \
                    (-Q2[z-1,t] + Q2[z+1,t])**2/(2.*dt**2) + \
                    (-Q3[z-1,t] + Q3[z+1,t])**2/(4.*dt**2) + \
                    (2*q0*(-Q2[z-1,t] + Q2[z+1,t])*Q4[z,t])/dt + \
                    (-Q4[z-1,t] + Q4[z+1,t])**2/(4.*dt**2) + \
                    4*q0*Q2[z,t]*((-Q1[z-1,t] + Q1[z+1,t])/(2.*dt) - \
                    (-Q4[z-1,t] + Q4[z+1,t])/(2.*dt)) + \
                    (2*q0*(-Q3[z-1,t] + Q3[z+1,t])*Q5[z,t])/dt - \
                    (2*q0*Q3[z,t]*(-Q5[z-1,t] + Q5[z+1,t]))/dt + \
                    (-Q5[z-1,t] + Q5[z+1,t])**2/(4.*dt**2)))/2.

                twisten += (kt*((-Q1[z-1,t] + Q1[z+1,t])**2/(4.*dt**2) - \
                    (2*q0*Q1[z,t]*(-Q2[z-1,t] + Q2[z+1,t]))/dt + \
                    (-Q2[z-1,t] + Q2[z+1,t])**2/(2.*dt**2) + \
                    (-Q3[z-1,t] + Q3[z+1,t])**2/(4.*dt**2) + \
                    (2*q0*(-Q2[z-1,t] + Q2[z+1,t])*Q4[z,t])/dt + \
                    (-Q4[z-1,t] + Q4[z+1,t])**2/(4.*dt**2) + \
                    4*q0*Q2[z,t]*((-Q1[z-1,t] + Q1[z+1,t])/(2.*dt) - \
                    (-Q4[z-1,t] + Q4[z+1,t])/(2.*dt)) + \
                    (2*q0*(-Q3[z-1,t] + Q3[z+1,t])*Q5[z,t])/dt - \
                    (2*q0*Q3[z,t]*(-Q5[z-1,t] + Q5[z+1,t]))/dt + \
                    (-Q5[z-1,t] + Q5[z+1,t])**2/(4.*dt**2)))/2.

                splay += (ks*((-Q3[z-1,t] + Q3[z+1,t])**2/(4.*dt**2) + \
                    (-(-Q1[z-1,t] + Q1[z+1,t])/(2.*dt) - (-Q4[z-1,t] + Q4[z+1,t])/(2.*dt))**2\
                    + (-Q5[z-1,t] + Q5[z+1,t])**2/(4.*dt**2)))/2.

                splayen += (ks*((-Q3[z-1,t] + Q3[z+1,t])**2/(4.*dt**2) + \
                    (-(-Q1[z-1,t] + Q1[z+1,t])/(2.*dt) - (-Q4[z-1,t] + Q4[z+1,t])/(2.*dt))**2\
                    + (-Q5[z-1,t] + Q5[z+1,t])**2/(4.*dt**2)))/2.

        enarrayinit.append(bulken + splayen + twisten + timeen + surfaceen)
        bulkenarr.append(bulken)
        twistenarr.append(twisten)
        splayenarr.append(splayen)
        timeenarr.append(timeen)
    #enarrayinit = np.array(enarrayinit)
    print(np.shape(enarrayinit))

    energy = (bulk + splay + twist + timeen + surface)# / (Lz*Lt)

    #calculategrad.calcgrad(guess,original,GradE,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk)
    np.savetxt("ks%ekt%eabc%e%e%esig%d/energyevolution.dat" % (ks,kt,a,b,c,sig), enarrayinit)
    np.savetxt("ks%ekt%eabc%e%e%esig%d/bulkevolution.dat" % (ks,kt,a,b,c,sig), bulkenarr)
    np.savetxt("ks%ekt%eabc%e%e%esig%d/twistevolution.dat" % (ks,kt,a,b,c,sig), twistenarr)
    np.savetxt("ks%ekt%eabc%e%e%esig%d/splayevolution.dat" % (ks,kt,a,b,c,sig), splayenarr)
    np.savetxt("ks%ekt%eabc%e%e%esig%d/timeevolution.dat" % (ks,kt,a,b,c,sig), timeenarr)

    return(energy)
