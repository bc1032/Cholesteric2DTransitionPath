#Libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize
import calculategrad

def calcenergy(guess,original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk):
    #print("Energy Called")
    dx = 1
    dz = 1
    dt = 1
    Qt1,Qt2,Qt3,Qt4,Qt5 = (2.0*s/3.0),0.0,0.0,-s/3.0,0.0
    #splay,twist,bend,surface,bulk = 0.0,0.0,0.0,0.0,0.0
    # timeen = 0.0
    enarrayfinal = []

    for t in range(0,Lt):
        for x in range(0,Lx):
            for z in range(0,Lz):
                Q1[x,z,t] = guess[5*x*Lz + 5*z + 5*t*Lx*Lz]
                Q2[x,z,t] = guess[5*x*Lz + 5*z + 5*t*Lx*Lz + 1]
                Q3[x,z,t] = guess[5*x*Lz + 5*z + 5*t*Lx*Lz + 2]
                Q4[x,z,t] = guess[5*x*Lz + 5*z + 5*t*Lx*Lz + 3]
                Q5[x,z,t] = guess[5*x*Lz + 5*z + 5*t*Lx*Lz + 4]
        enarrayinit = []

    for t in range(2,Lt-2):
        splayen,twisten,benden,surfaceen,bulken, timeenen = 0.0,0.0,0.0,0.0,0.0,0.0

        for x in range(0,Lx):

            if x == 0:
                xl = Lx-1
                xr = 1
            elif x == Lx-1:
                xl = Lx-2
                xr = 0
            else:
                xl = x - 1
                xr = x + 1
            for z in range(2,Lz-2):
                if t == 0:
                    Q1[x,z,t] = original[5*t*Lx + 5*z]
                    Q2[x,z,t] = original[5*t*Lx + 5*z + 1]
                    Q3[x,z,t] = original[5*t*Lx + 5*z + 2]
                    Q4[x,z,t] = original[5*t*Lx + 5*z + 3]
                    Q5[x,z,t] = original[5*t*Lx + 5*z + 4]
                if t == Lt-1:
                    Q1[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz]
                    Q2[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 1]
                    Q3[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 2]
                    Q4[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 3]
                    Q5[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 4]

                if t != 0 or t != Lt-1:
                    timeen += (sig*((Q1[x,z,t-1] - Q1[x,z,t+1])**2 + 2*(Q2[x,z,t-1] - Q2[x,z,t+1])**2\
                    + 2*(Q3[x,z,t-1] - Q3[x,z,t+1])**2 + (Q4[x,z,t-1] - Q4[x,z,t+1])**2 + \
                    (Q1[x,z,t-1] - Q1[x,z,t+1] + Q4[x,z,t-1] - Q4[x,z,t+1])**2 + 2*(Q5[x,z,t-1]\
                    - Q5[x,z,t+1])**2))/(8.*dt**2)
                #
                # if z == 0:
                #     Q1[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz]
                #     Q2[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 1]
                #     Q3[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 2]
                #     Q4[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 3]
                #     Q5[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 4]
                #
                #
                #     surface += 0.5*ws*((Q1[x,z,t]-Qt1)**2 + 2*(Q2[x,z,t]-Qt2)**2 \
                #     + 2*(Q3[x,z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[x,z,t] - Q4[x,z,t])**2 \
                #     + (Q4[x,z,t] - Qt4)**2 +2*(Q5[x,z,t] - Qt5)**2 )
                #
                #     # bulk += -(a*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                #     #         c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)**2 + \
                #     #         b*(Q1[x,z,t]**2*Q4[x,z,t] + (-Q2[x,z,t]**2 + Q3[x,z,t]**2)*Q4[x,z,t] - 2*Q2[x,z,t]*Q3[x,z,t]*Q5[x,z,t] +\
                #     #         Q1[x,z,t]*(-Q2[x,z,t]**2 + Q4[x,z,t]**2 + Q5[x,z,t]**2))
                # if z == Lz-1:
                #     Q1[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz]
                #     Q2[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 1]
                #     Q3[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 2]
                #     Q4[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 3]
                #     Q5[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 4]
                #
                #
                #     surface += 0.5*ws*( (Q1[x,z,t]-Qt1)**2 + 2*(Q2[x,z,t]-Qt2)**2 \
                #     + 2*(Q3[x,z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[x,z,t] - Q4[x,z,t])**2 \
                #     + (Q4[x,z,t] - Qt4)**2 +2*(Q5[x,z,t] - Qt5)**2 )
                #
                #     bulk += -(a*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                #             c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)**2 + \
                #             b*(Q1[x,z,t]**2*Q4[x,z,t] + (-Q2[x,z,t]**2 + Q3[x,z,t]**2)*Q4[x,z,t] - 2*Q2[x,z,t]*Q3[x,z,t]*Q5[x,z,t] +\
                #             Q1[x,z,t]*(-Q2[x,z,t]**2 + Q4[x,z,t]**2 + Q5[x,z,t]**2))
                else:

                    bulk += -(a*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                            c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)**2 + \
                            b*(Q1[x,z,t]**2*Q4[x,z,t] + (-Q2[x,z,t]**2 + Q3[x,z,t]**2)*Q4[x,z,t] - 2*Q2[x,z,t]*Q3[x,z,t]*Q5[x,z,t] +\
                            Q1[x,z,t]*(-Q2[x,z,t]**2 + Q4[x,z,t]**2 + Q5[x,z,t]**2))

                    twist += (kt*(dx**2*(4*dz*q0*Q1[x,z,t] + Q2[x,z01,t] - Q2[x,z+1,t])**2 + dz**2*(-Q2[xl,z,t] + Q2[xr,z,t] + 4*dx*q0*Q3[x,z,t])**2 +\
                            (dx*Q1[x,z-1,t] - dx*Q1[x,z+1,t] + dz*(-4*dx*q0*Q2[x,z,t] - Q3[xl,z,t] + Q3[xr,z,t]))**2 + \
                            dx**2*(4*dz*q0*Q2[x,z,t] + Q4[x,z-1,t] - Q4[x,z+1,t])**2 + dz**2*(-Q4[xl,z,t] + Q4[xr,z,t] + 4*dx*q0*Q5[x,z,t])**2 - \
                            32*dx**2*dz**2*q0**2*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2) + \
                            (dx*(-Q3[x,z-1,t] + Q3[x,z+1,t]) + dz*(-Q1[xl,z,t] + Q1[xr,z,t] - Q4[xl,z,t] + Q4[xr,z,t] + 4*dx*q0*Q5[x,z,t]))**2 + \
                            dx**2*(4*dz*q0*Q3[x,z,t] + Q5[x,z-1,t] - Q5[x,z+1,t])**2 + dz**2*(4*dx*q0*(Q1[x,z,t] + Q4[x,z,t]) + Q5[xl,z,t] - Q5[xr,z,t])**2 +\
                            (dx*Q2[x,z-1,t] - dx*Q2[x,z+1,t] + dz*(-4*dx*q0*Q4[x,z,t] - Q5[xl,z,t] + Q5[xr,z,t]))**2))/(8.*dx**2*dz**2)

                    splay += (ks*((dz*Q1[xl,z,t] - dz*Q1[xr,z,t] + dx*(Q3[x,z-1,t] - Q3[x,z+1,t]))**2 + \
                            (dz*(-Q3[xl,z,t] + Q3[xr,z,t]) + dx*(Q1[x,z-1,t] - Q1[x,z+1,t] + Q4[x,z-1,t] - Q4[x,z+1,t]))**2 +\
                            (dz*Q2[xl,z,t] - dz*Q2[xr,z,t] + dx*(Q5[x,z-1,t] - Q5[x,z+1,t]))**2))/(8.*dx**2*dz**2)

    energy = (bulk + splay + twist + timeen + surface) #/ (Lz*Lt)

    calculategrad.calcgrad(guess,original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk)

        #print(energy)
    return(energy,GradE)

def writeenergy(guess,original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk):
    #print("Energy Called")

    dz = 1
    dt = 1
    dx = 1
    Qt1,Qt2,Qt3,Qt4,Qt5= (2.0*s/3.0),0.0,0.0,-s/3.0,0.0
    # splay,twist,bend,surface,bulk = 0.0,0.0,0.0,0.0,0.0
    # timeen = 0.0
    enarrayfinal = []

    for t in range(0,Lt):
        for x in range(0,Lx):
            for z in range(0,Lz):
                Q1[x,z,t] = guess[5*x*Lz + 5*z + 5*t*Lx*Lz]
                Q2[x,z,t] = guess[5*x*Lz + 5*z + 5*t*Lx*Lz + 1]
                Q3[x,z,t] = guess[5*x*Lz + 5*z + 5*t*Lx*Lz + 2]
                Q4[x,z,t] = guess[5*x*Lz + 5*z + 5*t*Lx*Lz + 3]
                Q5[x,z,t] = guess[5*x*Lz + 5*z + 5*t*Lx*Lz + 4]
    enarrayinit = []
    bulkenarr, splayenarr, twistenarr, timeenarr = [],[],[],[]
    for t in range(0,Lt):
        splayen,twisten,benden,surfaceen,bulken, timeenen = 0.0,0.0,0.0,0.0,0.0,0.0
        for x in range(0,Lx):
            if x == 0:
                xl = Lx-1
                xr = 1
            elif x == Lx-1:
                xl = Lx-2
                xr = 0
            else:
                xl = x - 1
                xr = x + 1
            for z in range(0,Lz):
                if t == 0:
                    Q1[x,z,t] = original[5*x*Lz + 5*z]
                    Q2[x,z,t] = original[5*x*Lz + 5*z + 1]
                    Q3[x,z,t] = original[5*x*Lz + 5*z + 2]
                    Q4[x,z,t] = original[5*x*Lz + 5*z + 3]
                    Q5[x,z,t] = original[5*x*Lz + 5*z + 4]
                if t == Lt-1:
                    Q1[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz]
                    Q2[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 1]
                    Q3[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 2]
                    Q4[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 3]
                    Q5[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 4]

                else:
                    timeen += (sig*((-Q1[z,t-1] + Q1[z,t+1])**2/(4.*dz**2) + (-Q2[z,t-1] + Q2[z,t+1])**2/(2.*dz**2)\
                            + (-Q3[z,t-1] + Q3[z,t+1])**2/(2.*dz**2) + \
                            (-Q4[z,t-1] + Q4[z,t+1])**2/(4.*dz**2) + \
                            (-(-Q1[z,t-1] + Q1[z,t+1])/(2.*dz) - (-Q4[z,t-1] + Q4[z,t+1])/(2.*dz))**2 +\
                            (-Q5[z,t-1] + Q5[z,t+1])**2/(2.*dz**2)))/2.

                if z == 0:
                    Q1[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz]
                    Q2[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 1]
                    Q3[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 2]
                    Q4[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 3]
                    Q5[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 4]


                    surface += 0.5*ws*((Q1[x,z,t]-Qt1)**2 + 2*(Q2[x,z,t]-Qt2)**2 \
                    + 2*(Q3[x,z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[x,z,t] - Q4[x,z,t])**2 \
                    + (Q4[x,z,t] - Qt4)**2 +2*(Q5[x,z,t] - Qt5)**2 )

                    surfaceen += 0.5*ws*((Q1[x,z,t]-Qt1)**2 + 2*(Q2[x,z,t]-Qt2)**2 \
                    + 2*(Q3[x,z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[x,z,t] - Q4[x,z,t])**2 \
                    + (Q4[x,z,t] - Qt4)**2 +2*(Q5[x,z,t] - Qt5)**2 )

                    bulk += -(a*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2\
                            + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                             c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2\
                             + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)**2 + \
                             b*(Q1[x,z,t]**2*Q4[x,z,t] + (-Q2[x,z,t]**2 + Q3[x,z,t]**2)*Q4[x,z,t]\
                             - 2*Q2[x,z,t]*Q3[x,z,t]*Q5[x,z,t]\
                             + Q1[x,z,t]*(-Q2[x,z,t]**2 + Q4[x,z,t]**2 + Q5[x,z,t]**2))
                if z == Lz-1:
                    Q1[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz]
                    Q2[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 1]
                    Q3[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 2]
                    Q4[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 3]
                    Q5[x,z,t] = original[5*x*Lz + 5*z + 5*t*Lx*Lz + 4]


                    surface += 0.5*ws*( (Q1[x,z,t]-Qt1)**2 + 2*(Q2[x,z,t]-Qt2)**2 \
                    + 2*(Q3[x,z,t]-Qt3)**2+(Qt1 + Qt4 - Q1[x,z,t] - Q4[x,z,t])**2 \
                    + (Q4[x,z,t] - Qt4)**2 +2*(Q5[x,z,t] - Qt5)**2 )

                    bulk += -(a*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2\
                            + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                             c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2\
                             + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)**2 + \
                             b*(Q1[x,z,t]**2*Q4[x,z,t] + (-Q2[x,z,t]**2 + Q3[x,z,t]**2)*Q4[x,z,t]\
                             - 2*Q2[x,z,t]*Q3[x,z,t]*Q5[x,z,t]\
                             + Q1[x,z,t]*(-Q2[x,z,t]**2 + Q4[x,z,t]**2 + Q5[x,z,t]**2))

                else:

                    bulk += -(a*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2\
                            + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                             c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2\
                             + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)**2 + \
                             b*(Q1[x,z,t]**2*Q4[x,z,t] + (-Q2[x,z,t]**2 + Q3[x,z,t]**2)*Q4[x,z,t]\
                             - 2*Q2[x,z,t]*Q3[x,z,t]*Q5[x,z,t]\
                             + Q1[x,z,t]*(-Q2[x,z,t]**2 + Q4[x,z,t]**2 + Q5[x,z,t]**2))

                    twist += (kt*((-2*q0*Q1[x,z,t] + (-Q2[x,z-1,t] + Q2[x,z+1,t])/(2.*dz))**2\
                            + ((-Q2[xl,z,t] + Q2[xr,z,t])/(2.*dx) + 2*q0*Q3[x,z,t])**2 +
                            ((-Q1[x,z-1,t] + Q1[x,z+1,t])/(2.*dz) + 2*q0*Q2[x,z,t]\
                            - (-Q3[xl,z,t] + Q3[xr,z,t])/(2.*dx))**2 +
                            (-2*q0*Q2[x,z,t] + (-Q4[x,z-1,t] + Q4[x,z+1,t])/(2.*dz))**2\
                            + ((-Q4[xl,z,t] + Q4[xr,z,t])/(2.*dx) + 2*q0*Q5[x,z,t])**2 +
                            ((-Q1[xl,z,t] + Q1[xr,z,t])/(2.*dx) + (-Q3[x,z-1,t] + Q3[x,z+1,t])/(2.*dz)\
                            + (-Q4[xl,z,t] + Q4[xr,z,t])/(2.*dx) + 2*q0*Q5[x,z,t])**2 -\
                            8*q0**2*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t]\
                             + Q4[x,z,t]**2 + Q5[x,z,t]**2) +
                            (-2*q0*Q3[x,z,t] + (-Q5[x,z-1,t] + Q5[x,z+1,t])/(2.*dz))**2 + \
                            ((-Q2[x,z-1,t] + Q2[x,z+1,t])/(2.*dz) + 2*q0*Q4[x,z,t] - (-Q5[xl,z,t]\
                            + Q5[xr,z,t])/(2.*dx))**2 + \
                            (-2*q0*(Q1[x,z,t] + Q4[x,z,t]) + (-Q5[xl,z,t] + Q5[xr,z,t])/(2.*dx))**2))/2.

                    bulken += -(a*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2\
                            + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                             c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2\
                             + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)**2 + \
                             b*(Q1[x,z,t]**2*Q4[x,z,t] + (-Q2[x,z,t]**2 + Q3[x,z,t]**2)*Q4[x,z,t]\
                             - 2*Q2[x,z,t]*Q3[x,z,t]*Q5[x,z,t]\
                             + Q1[x,z,t]*(-Q2[x,z,t]**2 + Q4[x,z,t]**2 + Q5[x,z,t]**2))

                    twisten += (kt*((-2*q0*Q1[x,z,t] + (-Q2[x,z-1,t] + Q2[x,z+1,t])/(2.*dz))**2\
                            + ((-Q2[xl,z,t] + Q2[xr,z,t])/(2.*dx) + 2*q0*Q3[x,z,t])**2 +
                            ((-Q1[x,z-1,t] + Q1[x,z+1,t])/(2.*dz) + 2*q0*Q2[x,z,t]\
                            - (-Q3[xl,z,t] + Q3[xr,z,t])/(2.*dx))**2 +
                            (-2*q0*Q2[x,z,t] + (-Q4[x,z-1,t] + Q4[x,z+1,t])/(2.*dz))**2\
                            + ((-Q4[xl,z,t] + Q4[xr,z,t])/(2.*dx) + 2*q0*Q5[x,z,t])**2 +
                            ((-Q1[xl,z,t] + Q1[xr,z,t])/(2.*dx) + (-Q3[x,z-1,t] + Q3[x,z+1,t])/(2.*dz)\
                            + (-Q4[xl,z,t] + Q4[xr,z,t])/(2.*dx) + 2*q0*Q5[x,z,t])**2 -\
                            8*q0**2*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t]\
                             + Q4[x,z,t]**2 + Q5[x,z,t]**2) +
                            (-2*q0*Q3[x,z,t] + (-Q5[x,z-1,t] + Q5[x,z+1,t])/(2.*dz))**2 + \
                            ((-Q2[x,z-1,t] + Q2[x,z+1,t])/(2.*dz) + 2*q0*Q4[x,z,t] - (-Q5[xl,z,t]\
                            + Q5[xr,z,t])/(2.*dx))**2 + \
                            (-2*q0*(Q1[x,z,t] + Q4[x,z,t]) + (-Q5[xl,z,t] + Q5[xr,z,t])/(2.*dx))**2))/2.

                    splay += (ks*(((-Q1[xl,z,t] + Q1[xr,z,t])/(2.*dx) + (-Q3[x,z-1,t]\
                            + Q3[x,z+1,t])/(2.*dz))**2 + (-(-Q1[x,z-1,t] + Q1[x,z+1,t])/(2.*dz)\
                            + (-Q3[xl,z,t] + Q3[xr,z,t])/(2.*dx) - (-Q4[x,z-1,t]\
                            + Q4[x,z+1,t])/(2.*dz))**2 + ((-Q2[xl,z,t] + Q2[xr,z,t])/(2.*dx)\
                            + (-Q5[x,z-1,t] + Q5[x,z+1,t])/(2.*dz))**2))/2.

                    splayen += (ks*(((-Q1[xl,z,t] + Q1[xr,z,t])/(2.*dx) + (-Q3[x,z-1,t]\
                            + Q3[x,z+1,t])/(2.*dz))**2 + (-(-Q1[x,z-1,t] + Q1[x,z+1,t])/(2.*dz)\
                            + (-Q3[xl,z,t] + Q3[xr,z,t])/(2.*dx) - (-Q4[x,z-1,t]\
                            + Q4[x,z+1,t])/(2.*dz))**2 + ((-Q2[xl,z,t] + Q2[xr,z,t])/(2.*dx)\
                            + (-Q5[x,z-1,t] + Q5[x,z+1,t])/(2.*dz))**2))/2.

        enarrayinit.append(bulken + splayen + twisten + timeenen + surfaceen)
        bulkenarr.append(bulken)
        twistenarr.append(twisten)
        splayenarr.append(splayen)
        timeenarr.append(timeenen)
        #enarrayinit = np.array(enarrayinit)
        #print(np.shape(enarrayinit))

    energy = (bulk + splay + twist + timeen + surface)# / (Lz*Lt)

    calculategrad.calcgrad(guess,original,GradE,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk)
    np.savetxt("energyevolution.dat", enarrayinit)
    np.savetxt("bulkevolution.dat", bulkenarr)
    np.savetxt("twistevolution.dat", twistenarr)
    np.savetxt("splayevolution.dat", splayenarr)
    np.savetxt("timeevolution.dat" , timeenarr)

    return(energy,GradE)
