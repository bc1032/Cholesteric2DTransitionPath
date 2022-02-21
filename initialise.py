#Libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize

def initialise(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1, Q2, Q3, Q4, Q5, Lx, Length, Time, winding0, winding1, a, b, c):

    strength = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
    theta = math.pi/2.0

    for j in range(0,Time):
        for x in range(0,Lx):
            for i in range(0,Length):

                Q1[x,i,j] = strength*((-1.0/3.0)\
                        + (math.sin(theta)**2)*((Time-1-j)/(Time-1))*math.cos(winding0*i*math.pi/(Length-1))**2\
                        +   ((j)/(Time-1))*math.cos(winding1*i*math.pi/(Length-1))**2\
                        )

                Q2[x,i,j] = strength*(math.sin(theta)**2)*(\
                    ((Time-1-j)/(Time-1))*math.cos(winding0*i*math.pi/(Length-1))*math.sin(winding0*i*math.pi/(Length-1))\
                +   ((j)/(Time-1))*math.cos(winding1*i*math.pi/(Length-1))*math.sin(winding1*i*math.pi/(Length-1))\
                )

                Q3[x,i,j] = strength*(math.sin(theta)*math.cos(theta))*(\
                    ((Time-1-j)/(Time-1))*math.cos(winding0*i*math.pi/(Length-1))**2\
                +   ((j)/(Time-1))*math.cos(winding1*i*math.pi/(Length-1))**2\
                )

                Q4[x,i,j] = strength*((-1.0/3.0)\
                        + (math.sin(theta)**2)*((Time-1-j)/(Time-1))*math.sin(winding0*i*math.pi/(Length-1))**2\
                        +   ((j)/(Time-1))*math.sin(winding1*i*math.pi/(Length-1))**2\
                        )

                Q5[x,i,j] = strength*(math.sin(theta))*math.cos(theta)*(\
                    ((Time-1-j)/(Time-1))*math.sin(winding0*i*math.pi/(Length-1))\
                +   ((j)/(Time-1))*math.sin(winding1*i*math.pi/(Length-1))\
                )


    return(Q1, Q2, Q3, Q4, Q5)

def redefinerandom(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1, Q2, Q3, Q4, Q5, Lx, Length, Time, winding0, winding1, a, b, c):
    rana,ranb,ranc,rand,rane = [],[],[],[],[]
    #for j in range(5,Time-5):

    for x in range(0,Lx):
        n = 0.0

        for y in range(2,Length-2):
            zang = n*math.pi/(Length-4)

            rannum = np.random.uniform(-1.0,1.0)
            if rannum > 0:
                rana.append(np.random.normal(0.0,1.0)*math.sin(zang))
            else:
                rana.append(0.0)

            rannum = np.random.uniform(-1.0,1.0)
            if rannum > 0:
                ranb.append(np.random.normal(0.0,1.0)*math.sin(zang))
            else:
                ranb.append(0.0)
            rannum = np.random.uniform(-1.0,1.0)
            if rannum > 0:
                ranc.append(np.random.normal(0.0,1.0)*math.sin(zang))
            else:
                ranc.append(0.0)

            rannum = np.random.uniform(-1.0,1.0)
            if rannum > 0:
                rand.append(np.random.normal(0.0,1.0)*math.sin(zang))
            else:
                rand.append(0.0)

            rannum = np.random.uniform(-1.0,1.0)
            if rannum > 0:
                rane.append(np.random.normal(0.0,1.0)*math.sin(zang))
            else:
                rane.append(0.0)

            n+=1

    n = 0.0
    for j in range(5,Time-5):
        testang = n*math.pi/(Time-10)
        rancnt = 0
        for x in range(0,Lx):
            for i in range(2,Length-2):
                Q1[x,i,j] += math.sin(testang)*rana[rancnt]
                Q2[x,i,j] += math.sin(testang)*ranb[rancnt]
                Q3[x,i,j] += math.sin(testang)*ranc[rancnt]
                Q4[x,i,j] += math.sin(testang)*rand[rancnt]
                Q5[x,i,j] += math.sin(testang)*rane[rancnt]
                rancnt+=1


        n+=1
        #rancnt+=1
        #print(rancnt)

    return(Q1, Q2, Q3, Q4, Q5)

# def redefine2(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1,Q2,Q3,Q4,Q5, Length, Time, winding0, winding1, a, b, c):
#     Nt = timesplitf-timesplit
#     nt = 0
#     strength = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
#
#     for j in range(timesplit,timesplitf):
#     #if i >= lowz and i <= highz:
#         N = highz-lowz
#         n = lowz
#         theta = math.pi/2.0 + math.sin(nt*math.pi/((Nt-1)))*math.pi/2.0
#         transstate=np.loadtxt('uniaxial_q_V_alpha_eta_4.5_lambda_1_index_1.dat')
#         for i in range(0,Length):
#             Q1[x,i,j] = transstate[i*5]
#             Q2[x,i,j] transstate[i*5 + 1]
#             Q3[x,i,j] = transstate[i*5 + 2]
#             Q4[x,i,j] = transstate[i*5 + 3]
#             Q5[x,i,j] = transstate[i*5 + 4]
#             n += 1.0
#         nt += 1.0
#
#     return(Q1,Q2,Q3,Q4,Q5)

def redefine(splitt,splitz,lowz,highz,timesplit,timesplitf,Q1, Q2, Q3, Q4, Q5, Lx, Length, Time, winding0, winding1, a, b, c):

    Nt = timesplitf-timesplit
    nt = 0
    strength = (b + math.sqrt(b**2 + 24*a*c))/(4.0*c)
    for x in range(0,Lx):

        for j in range(timesplit,timesplitf):
        #if i >= lowz and i <= highz:
            N = highz-lowz
            n = lowz
            theta = math.pi/2.0 + math.sin(nt*math.pi/((Nt-1)))*math.pi/2.0

            for i in range(lowz,highz):
                Q1[x,i,j] = strength*((-1.0/3.0)\
                        + (math.sin(theta)**2)*((Time-1-j)/(Time-1))*math.cos(winding0*i*math.pi/(Length-1))**2\
                        +   ((j)/(Time-1))*math.cos(winding1*i*math.pi/(Length-1))**2\
                        )

                Q2[x,i,j] = strength*(math.sin(theta)**2)*(\
                    ((Time-1-j)/(Time-1))*math.cos(winding0*i*math.pi/(Length-1))*math.sin(winding0*i*math.pi/(Length-1))\
                +   ((j)/(Time-1))*math.cos(winding0*i*math.pi/(Length-1))*math.sin(winding1*i*math.pi/(Length-1))\
                )

                Q3[x,i,j] = strength*(math.sin(theta)*math.cos(theta))*(\
                    ((Time-1-j)/(Time-1))*math.cos(winding0*i*math.pi/(Length-1))**2\
                +   ((j)/(Time-1))*math.cos(winding1*i*math.pi/(Length-1))**2\
                )

                Q4[x,i,j] = strength*((-1.0/3.0)\
                        + (math.sin(theta)**2)*((Time-1-j)/(Time-1))*math.sin(winding0*i*math.pi/(Length-1))**2\
                        +   ((j)/(Time-1))*math.sin(winding1*i*math.pi/(Length-1))**2\
                        )

                Q5[x,i,j] = strength*(math.sin(theta))*math.cos(theta)*(\
                    ((Time-1-j)/(Time-1))*math.sin(winding0*i*math.pi/(Length-1))\
                +   ((j)/(Time-1))*math.sin(winding1*i*math.pi/(Length-1))\
                )
                n += 1.0
            nt += 1.0


    return(Q1, Q2, Q3, Q4, Q5)
