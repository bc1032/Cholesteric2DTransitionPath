import math
import numpy as np

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
                xr = 1
            elif x == Lx-1:
                xl = Lx-2
                xr = 0
            else:
                xl = x - 1
                xr = x + 1
            for z in range(2,Lz-2):
#Bulk Grad
                GradE[5*x*Lz + 5*z + 5*t*Lx*Lz] = (2*Q1[x,z,t] + Q4[x,z,t])*(-a + b*Q4[x,z,t] + 2*c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q3[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2)) + \
                                                    (b + 4*c*Q1[x,z,t] + 2*c*Q4[x,z,t])*Q5[x,z,t]**2 + (dx*kt*q0*Q2[x,z-1,t] - b*dx*dz*Q2[x,z,t]**2 - \
                                                    kt*q0*(dx*Q2[x,z+1,t] - dz*Q5[xl,z,t] + dz*Q5[xr,z,t]))/(dx*dz)

                GradE[5*x*Lz + 5*z + 5*t*Lx*Lz + 1] = 4*c*Q2[x,z,t]**3 + kt*q0*((Q3[xl,z,t] - Q3[xr,z,t])/dx + (-Q1[x,z-1,t] + Q1[x,z+1,t]\
                                                    + Q4[x,z-1,t] - Q4[x,z+1,t])/dz) - 2*b*Q3[x,z,t]*Q5[x,z,t] - \
                                                    2*Q2[x,z,t]*(a - 2*c*Q1[x,z,t]**2 + b*Q4[x,z,t] + Q1[x,z,t]*(b - 2*c*Q4[x,z,t])\
                                                    - 2*c*(Q3[x,z,t]**2 + Q4[x,z,t]**2 + Q5[x,z,t]**2))

                GradE[5*x*Lz + 5*z + 5*t*Lx*Lz + 2] = 4*c*Q3[x,z,t]**3 - 2*b*Q2[x,z,t]*Q5[x,z,t] + Q3[x,z,t]*(-2*a + 2*b*Q4[x,z,t] + \
                                                    4*c*(Q1[x,z,t]**2 + Q2[x,z,t]**2 + Q1[x,z,t]*Q4[x,z,t] + Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                                                    (kt*q0*(-(dz*Q2[xl,z,t]) + dz*Q2[xr,z,t] + dx*(Q5[x,z-1,t] - Q5[x,z+1,t])))/(dx*dz)

                GradE[5*x*Lz + 5*z + 5*t*Lx*Lz + 3] = 2*c*Q1[x,z,t]**3 + b*Q3[x,z,t]**2 - 2*a*Q4[x,z,t] + 4*c*Q3[x,z,t]**2*Q4[x,z,t]\
                                                    + 4*c*Q4[x,z,t]**3 - Q2[x,z,t]**2*(b - 4*c*Q4[x,z,t]) + Q1[x,z,t]**2*(b + 6*c*Q4[x,z,t])\
                                                    + 4*c*Q4[x,z,t]*Q5[x,z,t]**2 + Q1[x,z,t]*(-a + 2*b*Q4[x,z,t] + 2*c*(Q2[x,z,t]**2 + Q3[x,z,t]**2 + 3*Q4[x,z,t]**2 + Q5[x,z,t]**2)) + \
                                                    (kt*q0*(-(dx*Q2[x,z-1,t]) + dx*Q2[x,z+1,t] + 2*dz*Q5[xl,z,t] - 2*dz*Q5[xr,z,t]))/(dx*dz)

                GradE[5*x*Lz + 5*z + 5*t*Lx*Lz + 4] = (-(dz*kt*q0*Q1[xl,z,t]) + dz*kt*q0*Q1[xr,z,t] - dx*kt*q0*Q3[x,z-1,t] - 2*b*dx*dz*Q2[x,z,t]*Q3[x,z,t] + dx*kt*q0*Q3[x,z+1,t] - \
                                                    2*dz*kt*q0*Q4[xl,z,t] + 2*dz*kt*q0*Q4[xr,z,t] + 2*dx*dz*\
                                                    (-a + 2*c*Q1[x,z,t]**2 + Q1[x,z,t]*(b + 2*c*Q4[x,z,t]) + 2*c*(Q2[x,z,t]**2\
                                                    + Q3[x,z,t]**2 + Q4[x,z,t]**2))*Q5[x,z,t] + 4*c*dx*dz*Q5[x,z,t]**3)/(dx*dz)
    return(GradE)
