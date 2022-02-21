import math
import numpy as np
import scipy
import matplotlib.pyplot as pl
import matplotlib.colors as colors

def calcgrad(guess,original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk):
    dz = 1
    dt = 1
    dx = 1
    scale = 1.0#/((Lz*Lt))
    Qt1,Qt2,Qt3,Qt4,Qt5 = (2.0*s/3.0),0.0,0.0,-s/3.0,0.0

    #Calculate Bulk Gradient Energy Term and Elastic Energy Terms
    for t in range(2, Lt-2):
        for x in range(0,Lx):
            if x == 0:
                xl = Lx-1
                xr = x + 1
                xll = Lx-2
                xrr = x + 2
            elif x == Lx-1:
                xl = x-1
                xr = 0
                xll = x-2
                xrr = 1
            elif x == 1:
                xl = x-1
                xr = x+1
                xll = Lx-1
                xrr = 3
            elif x == Lx-2:
                xl = x-1
                xr = x+1
                xll = x-2
                xrr = 0
            else:
                xl = x - 1
                xr = x + 1
                xll = x - 2
                xrr = x + 2
            for z in range(2,Lz-2):
#Bulk Grad
                GradE[5*x*Lz + z + 5*t*Lx*Lz] = -(dt*dz*(ks + kt)*Q1[xll,z,t])/(4.*dx) - (dt*dx*(ks + kt)*Q1[x,z-2,t])/(4.*dz)\
                    - (dx*dz*sig*Q1[x,z,t-2])/(2.*dt) - 2*a*dt*dx*dz*Q1[x,z,t] + \
                    (dt*dx*ks*Q1[x,z,t])/(2.*dz) + (dt*dz*ks*Q1[x,z,t])/(2.*dx) + (dt*dx*kt*Q1[x,z,t])/(2.*dz) + (dt*dz*kt*Q1[x,z,t])/(2.*dx) + (dx*dz*sig*Q1[x,z,t])/dt + \
                    4*c*dt*dx*dz*Q1[x,z,t]**3 - (dx*dz*sig*Q1[x,z,t+2])/(2.*dt) - (dt*dx*ks*Q1[x,z+2,t])/(4.*dz) - (dt*dx*kt*Q1[x,z+2,t])/(4.*dz) - \
                    (dt*dz*ks*Q1[xrr,z,t])/(4.*dx) - (dt*dz*kt*Q1[xrr,z,t])/(4.*dx) + 2*dt*dx*kt*q0*Q2[x,z-1,t] - b*dt*dx*dz*Q2[x,z,t]**2\
                    + 4*c*dt*dx*dz*Q1[x,z,t]*Q2[x,z,t]**2 - \
                    2*dt*dx*kt*q0*Q2[x,z+1,t] + 4*c*dt*dx*dz*Q1[x,z,t]*Q3[x,z,t]**2 - (dt*dz*kt*Q4[xll,z,t])/(4.*dx) - (dt*dx*ks*Q4[x,z-2,t])/(4.*dz) - \
                    (dx*dz*sig*Q4[x,z,t-2])/(4.*dt) - a*dt*dx*dz*Q4[x,z,t] + (dt*dx*ks*Q4[x,z,t])/(2.*dz) + (dt*dz*kt*Q4[x,z,t])/(2.*dx)\
                    + (dx*dz*sig*Q4[x,z,t])/(2.*dt) + \
                    2*b*dt*dx*dz*Q1[x,z,t]*Q4[x,z,t] + 6*c*dt*dx*dz*Q1[x,z,t]**2*Q4[x,z,t] + 2*c*dt*dx*dz*Q2[x,z,t]**2*Q4[x,z,t] + 2*c*dt*dx*dz*Q3[x,z,t]**2*Q4[x,z,t] + \
                    b*dt*dx*dz*Q4[x,z,t]**2 + 6*c*dt*dx*dz*Q1[x,z,t]*Q4[x,z,t]**2 + 2*c*dt*dx*dz*Q4[x,z,t]**3 -\
                    (dx*dz*sig*Q4[x,z,t+2])/(4.*dt) - (dt*dx*ks*Q4[x,z+2,t])/(4.*dz) - \
                    (dt*dz*kt*Q4[xrr,z,t])/(4.*dx) + 2*dt*dz*kt*q0*Q5[xl,z,t] + b*dt*dx*dz*Q5[x,z,t]**2 + 4*c*dt*dx*dz*Q1[x,z,t]*Q5[x,z,t]**2 + \
                    2*c*dt*dx*dz*Q4[x,z,t]*Q5[x,z,t]**2 - 2*dt*dz*kt*q0*Q5[xr,z,t]

                GradE[5*x*Lz + z + 5*t*Lx*Lz + 1] = -2*dt*dx*kt*q0*Q1[x,z-1,t] + 2*dt*dx*kt*q0*Q1[x,z+1,t] - (dt*dz*ks*Q2[xll,z,t])/(4.*dx) - (dt*dz*kt*Q2[xll,z,t])/(4.*dx) -\
                    (dt*dx*kt*Q2[x,z-2,t])/(2.*dz) - (dx*dz*sig*Q2[x,z,t-2])/(2.*dt) - 2*a*dt*dx*dz*Q2[x,z,t] + (dt*dz*ks*Q2[x,z,t])/(2.*dx) + (dt*dx*kt*Q2[x,z,t])/dz +\
                    (dt*dz*kt*Q2[x,z,t])/(2.*dx) + (dx*dz*sig*Q2[x,z,t])/dt - 2*b*dt*dx*dz*Q1[x,z,t]*Q2[x,z,t] + 4*c*dt*dx*dz*Q1[x,z,t]**2*Q2[x,z,t] + 4*c*dt*dx*dz*Q2[x,z,t]**3 -\
                    (dx*dz*sig*Q2[x,z,t+2])/(2.*dt) - (dt*dx*kt*Q2[x,z+2,t])/(2.*dz) - (dt*dz*ks*Q2[xrr,z,t])/(4.*dx) - (dt*dz*kt*Q2[xrr,z,t])/(4.*dx) +\
                    2*dt*dz*kt*q0*Q3[xl,z,t] + 4*c*dt*dx*dz*Q2[x,z,t]*Q3[x,z,t]**2 - 2*dt*dz*kt*q0*Q3[xr,z,t] + 2*dt*dx*kt*q0*Q4[x,z-1,t] - 2*b*dt*dx*dz*Q2[x,z,t]*Q4[x,z,t] +\
                    4*c*dt*dx*dz*Q1[x,z,t]*Q2[x,z,t]*Q4[x,z,t] + 4*c*dt*dx*dz*Q2[x,z,t]*Q4[x,z,t]**2 - 2*dt*dx*kt*q0*Q4[x,z+1,t] - (dt*ks*Q5[xl,z-1,t])/4. +\
                    (dt*kt*Q5[xl,z-1,t])/4. + (dt*ks*Q5[xl,z+1,t])/4. - (dt*kt*Q5[xl,z+1,t])/4. - 2*b*dt*dx*dz*Q3[x,z,t]*Q5[x,z,t] +\
                    4*c*dt*dx*dz*Q2[x,z,t]*Q5[x,z,t]**2 + (dt*ks*Q5[xr,z-1,t])/4. - (dt*kt*Q5[xr,z-1,t])/4. + (dt*(-ks + kt)*Q5[xr,z+1,t])/4.

                GradE[5*x*Lz + z + 5*t*Lx*Lz + 2] = -2*dt*dz*kt*q0*Q2[xl,z,t] + 2*dt*dz*kt*q0*Q2[xr,z,t] - (dt*dz*ks*Q3[xll,z,t])/(4.*dx) - (dt*dz*kt*Q3[xll,z,t])/(4.*dx) -\
                    (dt*dx*ks*Q3[x,z-2,t])/(4.*dz) - (dt*dx*kt*Q3[x,z-2,t])/(4.*dz) - (dx*dz*sig*Q3[x,z,t-2])/(2.*dt) - 2*a*dt*dx*dz*Q3[x,z,t] + (dt*dx*ks*Q3[x,z,t])/(2.*dz) + \
                    (dt*dz*ks*Q3[x,z,t])/(2.*dx) + (dt*dx*kt*Q3[x,z,t])/(2.*dz) + (dt*dz*kt*Q3[x,z,t])/(2.*dx) + (dx*dz*sig*Q3[x,z,t])/dt + 4*c*dt*dx*dz*Q1[x,z,t]**2*Q3[x,z,t] + \
                    4*c*dt*dx*dz*Q2[x,z,t]**2*Q3[x,z,t] + 4*c*dt*dx*dz*Q3[x,z,t]**3 - (dx*dz*sig*Q3[x,z,t+2])/(2.*dt) - (dt*dx*ks*Q3[x,z+2,t])/(4.*dz) - \
                    (dt*dx*kt*Q3[x,z+2,t])/(4.*dz) - (dt*dz*ks*Q3[xrr,z,t])/(4.*dx) - (dt*dz*kt*Q3[xrr,z,t])/(4.*dx) + (dt*ks*Q4[xl,z-1,t])/4. - \
                    (dt*kt*Q4[xl,z-1,t])/4. - (dt*ks*Q4[xl,z+1,t])/4. + (dt*kt*Q4[xl,z+1,t])/4. + 2*b*dt*dx*dz*Q3[x,z,t]*Q4[x,z,t] + \
                    4*c*dt*dx*dz*Q1[x,z,t]*Q3[x,z,t]*Q4[x,z,t] + 4*c*dt*dx*dz*Q3[x,z,t]*Q4[x,z,t]**2 - (dt*ks*Q4[xr,z-1,t])/4. + (dt*kt*Q4[xr,z-1,t])/4. + \
                    (dt*ks*Q4[xr,z+1,t])/4. - (dt*kt*Q4[xr,z+1,t])/4. + 2*dt*dx*kt*q0*Q5[x,z-1,t] - 2*b*dt*dx*dz*Q2[x,z,t]*Q5[x,z,t] + \
                    4*c*dt*dx*dz*Q3[x,z,t]*Q5[x,z,t]**2 - 2*dt*dx*kt*q0*Q5[x,z+1,t]

                GradE[5*x*Lz + z + 5*t*Lx*Lz + 3] = -(dt*dz*kt*Q1[xll,z,t])/(4.*dx) - (dt*dx*ks*Q1[x,z-2,t])/(4.*dz) - \
                    (dx*dz*sig*Q1[x,z,t-2])/(4.*dt) - a*dt*dx*dz*Q1[x,z,t] + \
                    (dt*dx*ks*Q1[x,z,t])/(2.*dz) + (dt*dz*kt*Q1[x,z,t])/(2.*dx) + \
                    (dx*dz*sig*Q1[x,z,t])/(2.*dt) + b*dt*dx*dz*Q1[x,z,t]**2 + \
                    2*c*dt*dx*dz*Q1[x,z,t]**3 - (dx*dz*sig*Q1[x,z,t+2])/(4.*dt) -\
                    (dt*dx*ks*Q1[x,z+2,t])/(4.*dz) - (dt*dz*kt*Q1[xrr,z,t])/(4.*dx) -\
                    2*dt*dx*kt*q0*Q2[x,z-1,t] - b*dt*dx*dz*Q2[x,z,t]**2 + \
                    2*c*dt*dx*dz*Q1[x,z,t]*Q2[x,z,t]**2 + 2*dt*dx*kt*q0*Q2[x,z+1,t] +\
                    (dt*ks*Q3[xl,z-1,t])/4. - (dt*kt*Q3[xl,z-1,t])/4. - \
                    (dt*ks*Q3[xl,z+1,t])/4. + (dt*kt*Q3[xl,z+1,t])/4. + \
                    b*dt*dx*dz*Q3[x,z,t]**2 + 2*c*dt*dx*dz*Q1[x,z,t]*Q3[x,z,t]**2 - \
                    (dt*ks*Q3[xr,z-1,t])/4. + (dt*kt*Q3[xr,z-1,t])/4. + \
                    (dt*ks*Q3[xr,z+1,t])/4. - (dt*kt*Q3[xr,z+1,t])/4. - \
                    (dt*dz*kt*Q4[xll,z,t])/(2.*dx) - (dt*dx*ks*Q4[x,z-2,t])/(4.*dz) -\
                    (dt*dx*kt*Q4[x,z-2,t])/(4.*dz) - (dx*dz*sig*Q4[x,z,t-2])/(2.*dt) -\
                    2*a*dt*dx*dz*Q4[x,z,t] + (dt*dx*ks*Q4[x,z,t])/(2.*dz) + \
                    (dt*dx*kt*Q4[x,z,t])/(2.*dz) + (dt*dz*kt*Q4[x,z,t])/dx + \
                    (dx*dz*sig*Q4[x,z,t])/dt + 2*b*dt*dx*dz*Q1[x,z,t]*Q4[x,z,t] + \
                    6*c*dt*dx*dz*Q1[x,z,t]**2*Q4[x,z,t] + 4*c*dt*dx*dz*Q2[x,z,t]**2*Q4[x,z,t] + \
                    4*c*dt*dx*dz*Q3[x,z,t]**2*Q4[x,z,t] + 6*c*dt*dx*dz*Q1[x,z,t]*Q4[x,z,t]**2 + \
                    4*c*dt*dx*dz*Q4[x,z,t]**3 - (dx*dz*sig*Q4[x,z,t+2])/(2.*dt) - \
                    (dt*dx*ks*Q4[x,z+2,t])/(4.*dz) - (dt*dx*kt*Q4[x,z+2,t])/(4.*dz) - \
                    (dt*dz*kt*Q4[xrr,z,t])/(2.*dx) + 4*dt*dz*kt*q0*Q5[xl,z,t] + \
                    2*c*dt*dx*dz*Q1[x,z,t]*Q5[x,z,t]**2 + 4*c*dt*dx*dz*Q4[x,z,t]*Q5[x,z,t]**2 - \
                    4*dt*dz*kt*q0*Q5[xr,z,t]

                GradE[5*x*Lz + z + 5*t*Lx*Lz + 4] =         -2*dt*dz*kt*q0*Q1[xl,z,t] + 2*dt*dz*kt*q0*Q1[xr,z,t] - \
                    (dt*ks*Q2[xl,z-1,t])/4. + (dt*kt*Q2[xl,z-1,t])/4. +\
                    (dt*ks*Q2[xl,z+1,t])/4. - (dt*kt*Q2[xl,z+1,t])/4. + \
                    (dt*ks*Q2[xr,z-1,t])/4. - (dt*kt*Q2[xr,z-1,t])/4. - \
                    (dt*ks*Q2[xr,z+1,t])/4. + (dt*kt*Q2[xr,z+1,t])/4. - \
                    2*dt*dx*kt*q0*Q3[x,z-1,t] - 2*b*dt*dx*dz*Q2[x,z,t]*Q3[x,z,t] + \
                    2*dt*dx*kt*q0*Q3[x,z+1,t] - 4*dt*dz*kt*q0*Q4[xl,z,t] + \
                    4*dt*dz*kt*q0*Q4[xr,z,t] - (dt*dz*kt*Q5[xll,z,t])/(2.*dx) - \
                    (dt*dx*ks*Q5[x,z-2,t])/(4.*dz) - (dt*dx*kt*Q5[x,z-2,t])/(4.*dz) - \
                    (dx*dz*sig*Q5[x,z,t-2])/(2.*dt) - 2*a*dt*dx*dz*Q5[x,z,t] + \
                    (dt*dx*ks*Q5[x,z,t])/(2.*dz) + (dt*dx*kt*Q5[x,z,t])/(2.*dz) + \
                    (dt*dz*kt*Q5[x,z,t])/dx + (dx*dz*sig*Q5[x,z,t])/dt + \
                    2*b*dt*dx*dz*Q1[x,z,t]*Q5[x,z,t] + 4*c*dt*dx*dz*Q1[x,z,t]**2*Q5[x,z,t] + \
                    4*c*dt*dx*dz*Q2[x,z,t]**2*Q5[x,z,t] + 4*c*dt*dx*dz*Q3[x,z,t]**2*Q5[x,z,t] + \
                    4*c*dt*dx*dz*Q1[x,z,t]*Q4[x,z,t]*Q5[x,z,t] + \
                    4*c*dt*dx*dz*Q4[x,z,t]**2*Q5[x,z,t] + 4*c*dt*dx*dz*Q5[x,z,t]**3 - \
                    (dx*dz*sig*Q5[x,z,t+2])/(2.*dt) - (dt*dx*ks*Q5[x,z+2,t])/(4.*dz) - \
                    (dt*dx*kt*Q5[x,z+2,t])/(4.*dz) - (dt*dz*kt*Q5[xrr,z,t])/(2.*dx)
    return(GradE)

def calcgradold(guess,original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk):
    dz = 1
    dt = 1
    dx = 1
    scale = 1.0#/((Lz*Lt))
    Qt1,Qt2,Qt3,Qt4,Qt5 = (2.0*s/3.0),0.0,0.0,-s/3.0,0.0

    #Calculate Bulk Gradient Energy Term and Elastic Energy Terms
    for t in range(2, Lt-2):
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
#Bulk Grad
                GradE[5*x*Lz + z + 5*t*Lx*Lz] = (2*Q1[x,z,t] + Q4[x,z,t])*(-a + b*Q4[x,z,t] + 2*c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2)) + \
                                                    (b + 4*c*Q1[x,z,t] + 2*c*Q4[x,z,t])*Q5[x,z,t]**2 + (dx*kt*q0*Q2[x,z-1,t] - b*dx*dz*Q2[x,z,t]**2 - \
                                                    kt*q0*(dx*Q2[x,z+1,t] - dz*Q5[xl,z,t] + dz*Q5[xr,z,t]))/(dx*dz)

                GradE[5*x*Lz + z + 5*t*Lx*Lz + 1] = 4*c*Q2[x,z,t]**3 + kt*q0*((Q3[xl,z,t] - Q3[xr,z,t])/dx + (-Q1[x,z-1,t] + Q1[x,z+1,t]\
                                                    + Q4[x,z-1,t] - Q4[x,z+1,t])/dz) - 2*b*Q3[x,z,t]*Q5[x,z,t] - \
                                                    2*Q2[x,z,t]*(a - 2*c*Q1[x,z,t]**2 + b*Q4[x,z,t] + Q1[x,z,t]*(b - 2*c*Q4[x,z,t])\
                                                    - 2*c*(Q3[x,z,t]**2 + Q4[x,z,t]**2 + Q5[x,z,t]**2))

                GradE[5*x*Lz + z + 5*t*Lx*Lz + 2] = 4*c*Q3[x,z,t]**3 - 2*b*Q2[x,z,t]*Q5[x,z,t] + Q3[x,z,t]*(-2*a + 2*b*Q4[x,z,t] + \
                                                    4*c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                                                    (kt*q0*(-(dz*Q2[xl,z,t]) + dz*Q2[xr,z,t] + dx*(Q5[x,z-1,t] - Q5[x,z+1,t])))/(dx*dz)

                GradE[5*x*Lz + z + 5*t*Lx*Lz + 3] = 2*c*Q1[x,z,t]**3 + b*Q3[x,z,t]**2 - 2*a*Q4[x,z,t] + 4*c*Q3[x,z,t]**2*Q4[x,z,t]\
                                                    + 4*c*Q4[x,z,t]**3 - Q2[x,z,t]**2*(b - 4*c*Q4[x,z,t]) + Q1[x,z,t]**2*(b + 6*c*Q4[x,z,t])\
                                                    + 4*c*Q4[x,z,t]*Q5[x,z,t]**2 + Q1[x,z,t]*(-a + 2*b*Q4[x,z,t] + 2*c*(Q2[x,z,t]**2 + Q3[x,z,t]**2 + 3*Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                                                    (kt*q0*(-(dx*Q2[x,z-1,t]) + dx*Q2[x,z+1,t] + 2*dz*Q5[xl,z,t] - 2*dz*Q5[xr,z,t]))/(dx*dz)

                GradE[5*x*Lz + z + 5*t*Lx*Lz + 4] = (-(dz*kt*q0*Q1[xl,z,t]) + dz*kt*q0*Q1[xr,z,t] - dx*kt*q0*Q3[x,z-1,t] - 2*b*dx*dz*Q2[x,z,t]*Q3[x,z,t] + dx*kt*q0*Q3[x,z+1,t] - \
                                                    2*dz*kt*q0*Q4[xl,z,t] + 2*dz*kt*q0*Q4[xr,z,t] + 2*dx*dz*\
                                                    (-a + 2*c*Q1[x,z,t]**2 + Q1[x,z,t]*(b + 2*c*Q4[x,z,t]) + 2*c*(Q2[x,z,t]**2\
                                                    + Q3[x,z,t]**2 + Q4[x,z,t]**2))*Q5[x,z,t] + 4*c*dx*dz*Q5[x,z,t]**3)/(dx*dz)
    return(GradE)

def Hessian(sig,Lx,Lz,Lt,ks,kt,q0,s,a,b,c,Q1,Q2,Q3,Q4,Q5):
    dz = 1
    dt = 1
    dx = 1
    scale = 1.0#/((Lz*Lt))
    Qt1,Qt2,Qt3,Qt4,Qt5 = (2.0*s/3.0),0.0,0.0,-s/3.0,0.0
    Hessian = np.zeros((5*Lx*Lz,Lx*Lz*5))
    Hessian2=[]
    t = 2
    #Calculate Bulk Gradient Energy Term and Elastic Energy Terms

    for param in range(0,(Lx*Lz)):
        #if param % 200 == 0:
        #    print(param)
        for onsitex in range(0,Lx):
            for onsitez in range(0,Lz):
                for x in range(0,Lx):
                    onsitez = 0
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
                        zb = z - 1
                        za = z + 1
                        #if z % 10 == 0:
                            #print(z)
                            #test = 0
                        if onsitex == x and onsitez == z:
                            rightQ1Q5 = -9*dx*dz**2*kt*q0
                            rightQ2Q3 = -9*dx*dz**2*kt*q0
                            rightQ3Q2 = 9*dx*dz**2*kt*q0
                            rightQ4Q5 = -18*dx*dz**2*kt*q0
                            rightQ5Q1 = 9*dx*dz**2*kt*q0
                            rightQ5Q4 = 18*dx*dz**2*kt*q0

                            leftQ1Q5 = 9*dx*dz**2*kt*q0
                            leftQ2Q3 = 9*dx*dz**2*kt*q0
                            leftQ3Q2 = -9*dx*dz**2*kt*q0
                            leftQ4Q5 = 18*dx*dz**2*kt*q0
                            leftQ5Q1 = -9*dx*dz**2*kt*q0
                            leftQ5Q4 = -18*dx*dz**2*kt*q0

                            onsiteQ1Q1 = 18*dx**2*dz**2*(-a + b*Q4[x,z,t] + 3*c*(2*Q1[x,z,t]**2\
                                        + 2*Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2) + 2*c*(Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q5[x,z,t]**2))
                            onsiteQ1Q2 = -18*dx**2*dz**2*Q2[x,z,t]*(b - 2*c*(2*Q1[x,z,t] + Q4[x,z,t]))
                            onsiteQ1Q3 = 36*c*dx**2*dz**2*Q3[x,z,t]*(2*Q1[x,z,t] + Q4[x,z,t])
                            onsiteQ1Q4 = 9*dx**2*dz**2*(-a + 6*c*Q1[x,z,t]**2 + 2*b*Q4[x,z,t]\
                                         + 2*Q1[x,z,t]*(b + 6*c*Q4[x,z,t]) + 2*c*(Q2[x,z,t]**2 + Q3[x,z,t]**2 + 3*Q4[x,z,t]**2 + Q5[x,z,t]**2))
                            onsiteQ1Q5 = 18*dx**2*dz**2*(b + 4*c*Q1[x,z,t] + 2*c*Q4[x,z,t])*Q5[x,z,t]
                            onsiteQ2Q2 = 9*dx**2*dz**2*(4*c*Q1[x,z,t]**2 - 2*(a + b*Q4[x,z,t])\
                                        - 2*Q1[x,z,t]*(b - 2*c*Q4[x,z,t]) + 4*c*(3*Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q4[x,z,t]**2 + Q5[x,z,t]**2))
                            onsiteQ2Q3 = 9*dx**2*dz**2*(8*c*Q2[x,z,t]*Q3[x,z,t] - 2*b*Q5[x,z,t])
                            onsiteQ2Q4 = -18*dx**2*dz**2*Q2[x,z,t]*(b - 2*c*(Q1[x,z,t] + 2*Q4[x,z,t]))
                            onsiteQ2Q5 = 9*dx**2*dz**2*(-2*b*Q3[x,z,t] + 8*c*Q2[x,z,t]*Q5[x,z,t])
                            onsiteQ3Q3 = 18*dx**2*dz**2*(-a + b*Q4[x,z,t] + 2*c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + 3*Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2))
                            onsiteQ3Q4 = 18*dx**2*dz**2*Q3[x,z,t]*(b + 2*c*Q1[x,z,t] + 4*c*Q4[x,z,t])
                            onsiteQ3Q5 = -18*dx**2*dz**2*(b*Q2[x,z,t] - 4*c*Q3[x,z,t]*Q5[x,z,t])
                            onsiteQ4Q4 = 9*dx**2*dz**2*(-2*a + 6*c*Q1[x,z,t]**2 + 2*Q1[x,z,t]*(b + 6*c*Q4[x,z,t]) + 4*c*(Q2[x,z,t]**2 + Q3[x,z,t]**2 + 3*Q4[x,z,t]**2 + Q5[x,z,t]**2))
                            onsiteQ4Q5 = 36*c*dx**2*dz**2*(Q1[x,z,t] + 2*Q4[x,z,t])*Q5[x,z,t]
                            onsiteQ5Q5 = 18*dx**2*dz**2*(-a + 2*c*Q1[x,z,t]**2 + Q1[x,z,t]*(b + 2*c*Q4[x,z,t]) + 2*c*(Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q4[x,z,t]**2 + 3*Q5[x,z,t]**2))




                            Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*z)  ] =onsiteQ1Q1
                            Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*z) + 1  ] =onsiteQ1Q2
                            Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*z) + 2  ] =onsiteQ1Q3
                            Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*z) + 3  ] =onsiteQ1Q4
                            Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*z) + 4  ] =onsiteQ1Q5
                            Hessian[(5*x*Lz + 5*z)+1,(5*x*Lz + 5*z  )  ] =onsiteQ1Q2
                            Hessian[(5*x*Lz + 5*z)+1,(5*x*Lz + 5*z  ) + 1] =onsiteQ2Q2
                            Hessian[(5*x*Lz + 5*z)+1,(5*x*Lz + 5*z  ) + 2] =onsiteQ2Q3
                            Hessian[(5*x*Lz + 5*z)+1,(5*x*Lz + 5*z  ) + 3] =onsiteQ2Q4
                            Hessian[(5*x*Lz + 5*z)+1,(5*x*Lz + 5*z  ) + 4] =onsiteQ2Q5
                            Hessian[(5*x*Lz + 5*z)+2,(5*x*Lz + 5*z  )] =onsiteQ1Q3
                            Hessian[(5*x*Lz + 5*z)+2,(5*x*Lz + 5*z  ) + 1] =onsiteQ2Q3
                            Hessian[(5*x*Lz + 5*z)+2,(5*x*Lz + 5*z  ) + 2] =onsiteQ3Q3
                            Hessian[(5*x*Lz + 5*z)+2,(5*x*Lz + 5*z  ) + 3] =onsiteQ3Q4
                            Hessian[(5*x*Lz + 5*z)+2,(5*x*Lz + 5*z  ) + 4] =onsiteQ3Q5
                            Hessian[(5*x*Lz + 5*z)+3,(5*x*Lz + 5*z  ) + 0] =onsiteQ1Q4
                            Hessian[(5*x*Lz + 5*z)+3,(5*x*Lz + 5*z  ) + 1] =onsiteQ2Q4
                            Hessian[(5*x*Lz + 5*z)+3,(5*x*Lz + 5*z  ) + 2] =onsiteQ3Q4
                            Hessian[(5*x*Lz + 5*z)+3,(5*x*Lz + 5*z  ) + 3] =onsiteQ4Q4
                            Hessian[(5*x*Lz + 5*z)+3,(5*x*Lz + 5*z  ) + 4] =onsiteQ4Q5
                            Hessian[(5*x*Lz + 5*z)+4,(5*x*Lz + 5*z  ) + 0] =onsiteQ1Q5
                            Hessian[(5*x*Lz + 5*z)+4,(5*x*Lz + 5*z  ) + 1] =onsiteQ2Q5
                            Hessian[(5*x*Lz + 5*z)+4,(5*x*Lz + 5*z  ) + 2] =onsiteQ3Q5
                            Hessian[(5*x*Lz + 5*z)+4,(5*x*Lz + 5*z  ) + 3] =onsiteQ4Q5
                            Hessian[(5*x*Lz + 5*z)+4,(5*x*Lz + 5*z  ) + 4] =onsiteQ5Q5

                            Hessian[(5*x*Lz + 5*z),(5*xr*Lz + 5*z ) + 4] =rightQ1Q5
                            Hessian[(5*x*Lz + 5*z)+1,(5*xr*Lz + 5*z ) + 2] =rightQ2Q3
                            Hessian[(5*x*Lz + 5*z)+2,(5*xr*Lz + 5*z ) + 1] =rightQ3Q2
                            Hessian[(5*x*Lz + 5*z)+3,(5*xr*Lz + 5*z ) + 4] =rightQ4Q5
                            Hessian[(5*x*Lz + 5*z)+4,(5*xr*Lz + 5*z ) + 0] =rightQ5Q1
                            Hessian[(5*x*Lz + 5*z)+4,(5*xr*Lz + 5*z ) + 3] =rightQ5Q4

                            Hessian[(5*x*Lz + 5*z),(5*(xl)*Lz + 5*z ) + 4] =leftQ1Q5
                            Hessian[(5*x*Lz + 5*z)+1,(5*(xl)*Lz + 5*z ) + 2] =leftQ2Q3
                            Hessian[(5*x*Lz + 5*z)+2,(5*(xl)*Lz + 5*z ) + 1] =leftQ3Q2
                            Hessian[(5*x*Lz + 5*z)+3,(5*(xl)*Lz + 5*z ) + 4] =leftQ4Q5
                            Hessian[(5*x*Lz + 5*z)+4,(5*(xl)*Lz + 5*z ) + 0] =leftQ5Q1
                            Hessian[(5*x*Lz + 5*z)+4,(5*(xl)*Lz + 5*z ) + 3] =leftQ5Q4

                            if z!= Lz-1:
                                aboveQ1Q2 = -9*dx**2*dz*kt*q0
                                aboveQ2Q1 = 9*dx**2*dz*kt*q0
                                aboveQ2Q4 = -9*dx**2*dz*kt*q0
                                aboveQ3Q5 = -9*dx**2*dz*kt*q0
                                aboveQ4Q2 = 9*dx**2*dz*kt*q0
                                aboveQ5Q3 = 9*dx**2*dz*kt*q0

                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*za) + 1] = aboveQ1Q2
                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*za)] = aboveQ2Q1
                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*za) + 3] = aboveQ2Q4
                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*za) + 4] = aboveQ3Q5
                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*za) + 1] = aboveQ4Q2
                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*za) + 2] = aboveQ5Q3
                            if z!= 0:
                                belowQ1Q2 = 9*dx**2*dz*kt*q0
                                belowQ2Q1 = -9*dx**2*dz*kt*q0
                                belowQ2Q4 = 9*dx**2*dz*kt*q0
                                belowQ3Q5 = 9*dx**2*dz*kt*q0
                                belowQ4Q2 = -9*dx**2*dz*kt*q0
                                belowQ5Q3 = -9*dx**2*dz*kt*q0

                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*zb) + 1] = belowQ1Q2
                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*zb)] = belowQ2Q1
                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*zb) + 3] = belowQ2Q4
                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*zb) + 4] = belowQ3Q5
                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*zb) + 1] = belowQ4Q2
                                Hessian[(5*x*Lz + 5*z),(5*x*Lz + 5*zb) + 2] = belowQ5Q3

    pl.imshow(Hessian)
    pl.savefig('hessian.pdf')
    #pl.show()
    pl.close()
    np.savetxt("Hessian.dat",(Hessian[:,0]))

    a = np.asarray(np.linalg.eigvals(Hessian))
    #print(a)
    evalp, evaln = 0,0
    for i in range(0,len(a)):
        if a[i] > 0:
            evalp += 1
        elif a[i] < 0:
            evaln += 1

    np.savetxt("HessianEvals.dat", a)
    # pl.plot(len,evalp, label = 'Positive')
    # pl.plot(len,evaln, label = 'Negative')
    # pl.legend(loc="upper right")
    # pl.xlabel('$k_s, k_t = 1^{x}$')
    # pl.ylabel('Number of eigenvalues of Hessian')
    # pl.show()
    # pl.savefig('hessianevals.pdf')
    #return(Hessian)
    return(evalp,evaln)
