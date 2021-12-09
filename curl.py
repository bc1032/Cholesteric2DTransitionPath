import numpy as np
from numpy import random
import pylab as pl
import os
def curl(lz,lx,lt):
    lz, lx = 30,30
        lt = 200
    M = []
    data = np.loadtxt('energyarray.dat')
    #data = np.loadtxt('ftest.csv')
    #data = np.(lz,lx)
    i=1
    n=0

    if not os.path.exists('curl'):
        os.makedirs('curl')

    if not os.path.exists('curlvids'):
        os.makedirs('curlvids')

    if not os.path.exists('curlgraphs'):
        os.makedirs('curlgraphs')

    fcurl = open("curl/curl.dat",'w')

    for t in range(0,lt):
        if t % 20 == 0:
            print(t)
        for i in range(0,lx):
            for k in range(0,lz):

                if i == 0:
                    xl = lx-1
                    xr = i + 1
                elif i == lx-1:
                    xl = lx-2
                    xr = i-1
                else:
                    xl = i-1
                    xr = i+1

                za = k+1
                zb = k-1

                indexz = np.mod(int(i), lz)
                indexx = int(int(i)/lz)
                #print(5*i*lx + 5*k + 5*t*lx*lz)
                Q1 = data[5*i*lz + 5*k + 5*t*lx*lz]
                Q2 = data[5*i*lz + 5*k + 5*t*lx*lz + 1]
                Q3 = data[5*i*lz + 5*k + 5*t*lx*lz + 2]
                Q4 = data[5*i*lz + 5*k + 5*t*lx*lz + 3]
                Q5 = data[5*i*lz + 5*k + 5*t*lx*lz + 4]
                curlq = 0.0
                M = np.array([[Q1,Q2,Q3],[Q2,Q4,Q5],[Q3,Q5,-Q1-Q4]])

                #Q1R = data[5*xr*lz + 5*k + 5*t*lx*lz]

                Q1R = data[5*xr*lz + 5*k + 5*t*lx*lz]
                Q2R = data[5*xr*lz + 5*k + 5*t*lx*lz + 1]
                Q3R = data[5*xr*lz + 5*k + 5*t*lx*lz + 2]
                Q4R = data[5*xr*lz + 5*k + 5*t*lx*lz + 3]
                Q5R = data[5*xr*lz + 5*k + 5*t*lx*lz + 4]

                Q1L = data[5*xl*lz + 5*k + 5*t*lx*lz]
                Q2L = data[5*xl*lz + 5*k + 5*t*lx*lz + 1]
                Q3L = data[5*xl*lz + 5*k + 5*t*lx*lz + 2]
                Q4L = data[5*xl*lz + 5*k + 5*t*lx*lz + 3]
                Q5L = data[5*xl*lz + 5*k + 5*t*lx*lz + 4]

                ML = np.array([[Q1L,Q2L,Q3L],[Q2L,Q4L,Q5L],[Q3L,Q5L,-Q1L-Q4L]])
                MR = np.array([[Q1R,Q2R,Q3R],[Q2R,Q4R,Q5R],[Q3R,Q5R,-Q1R-Q4R]])

                lciv = np.array([[[0,0,0],[0,0,1],[0,-1,0]], [[0,0,-1],[0,0,0],[1,0,0]],[[0,1,0],[-1,0,0],[0,0,0]]])
                if k == 0:
                    Q1F = -(3.0/2.0)*data[5*i*lz + 5*k + 5*t*lx*lz] + 2.0*data[5*i*lz + 5*(k+1) + 5*t*lx*lz] - (1.0/2.0)*data[5*i*lz + 5*(k+2) + 5*t*lx*lz]
                    Q2F = -(3.0/2.0)*data[5*i*lz + 5*k + 5*t*lx*lz + 1] + 2.0*data[5*i*lz + 5*(k+1) + 5*t*lx*lz + 1] - (1.0/2.0)*data[5*i*lz + 5*(k+2) + 5*t*lx*lz + 1]
                    Q3F = -(3.0/2.0)*data[5*i*lz + 5*k + 5*t*lx*lz + 2] + 2.0*data[5*i*lz + 5*(k+1) + 5*t*lx*lz + 2] - (1.0/2.0)*data[5*i*lz + 5*(k+2) + 5*t*lx*lz + 2]
                    Q4F = -(3.0/2.0)*data[5*i*lz + 5*k + 5*t*lx*lz + 3] + 2.0*data[5*i*lz + 5*(k+1) + 5*t*lx*lz + 3] - (1.0/2.0)*data[5*i*lz + 5*(k+2) + 5*t*lx*lz + 3]
                    Q5F = -(3.0/2.0)*data[5*i*lz + 5*k + 5*t*lx*lz + 4] + 2.0*data[5*i*lz + 5*(k+1) + 5*t*lx*lz + 4] - (1.0/2.0)*data[5*i*lz + 5*(k+2) + 5*t*lx*lz + 4]
                    MZ = np.array([[Q1F,Q2F,Q3F],[Q2F,Q4F,Q5F],[Q3F,Q5F,-Q1F-Q4F]])
                elif k == lz - 1:
                    Q1B = (3.0/2.0)*data[5*i*lz + 5*k + 5*t*lx*lz] - 2.0*data[5*i*lz + 5*(k-1) + 5*t*lx*lz] + (1.0/2.0)*data[5*i*lz + 5*(k-2) + 5*t*lx*lz]
                    Q2B = (3.0/2.0)*data[5*i*lz + 5*k + 5*t*lx*lz + 1] - 2.0*data[5*i*lz + 5*(k-1) + 5*t*lx*lz + 1] + (1.0/2.0)*data[5*i*lz + 5*(k-2) + 5*t*lx*lz + 1]
                    Q3B = (3.0/2.0)*data[5*i*lz + 5*k + 5*t*lx*lz + 2] - 2.0*data[5*i*lz + 5*(k-1) + 5*t*lx*lz + 2] + (1.0/2.0)*data[5*i*lz + 5*(k-2) + 5*t*lx*lz + 2]
                    Q4B = (3.0/2.0)*data[5*i*lz + 5*k + 5*t*lx*lz + 3] - 2.0*data[5*i*lz + 5*(k-1) + 5*t*lx*lz + 3] + (1.0/2.0)*data[5*i*lz + 5*(k-2) + 5*t*lx*lz + 3]
                    Q5B = (3.0/2.0)*data[5*i*lz + 5*k + 5*t*lx*lz + 4] - 2.0*data[5*i*lz + 5*(k-1) + 5*t*lx*lz + 4] + (1.0/2.0)*data[5*i*lz + 5*(k-2) + 5*t*lx*lz + 4]
                    MZ = np.array([[Q1B,Q2B,Q3B],[Q2B,Q4B,Q5B],[Q3B,Q5B,-Q1B-Q4B]])
                else:
                    Q1A = data[5*i*lz + 5*za + 5*t*lx*lz]
                    Q2A = data[5*i*lz + 5*za + 5*t*lx*lz + 1]
                    Q3A = data[5*i*lz + 5*za + 5*t*lx*lz + 2]
                    Q4A = data[5*i*lz + 5*za + 5*t*lx*lz + 3]
                    Q5A = data[5*i*lz + 5*za + 5*t*lx*lz + 4]

                    Q1B = data[5*i*lz + 5*zb + 5*t*lx*lz]
                    Q2B = data[5*i*lz + 5*zb + 5*t*lx*lz + 1]
                    Q3B = data[5*i*lz + 5*zb + 5*t*lx*lz + 2]
                    Q4B = data[5*i*lz + 5*zb + 5*t*lx*lz + 3]
                    Q5B = data[5*i*lz + 5*zb + 5*t*lx*lz + 4]

                    MA = np.array([[Q1A,Q2A,Q3A],[Q2A,Q4A,Q5A],[Q3A,Q5A,-Q1A-Q4A]])
                    MB = np.array([[Q1B,Q2B,Q3B],[Q2B,Q4B,Q5B],[Q3B,Q5B,-Q1B-Q4B]])
                    MZ = np.subtract(MA, MB)
                MX = np.subtract(MR, ML)
                #if t == 0 or t == lt-1:
                curlq = 0.0
                #else:
                for qi in range(0,3):
                    for qj in range(0,3):
                        for ql in range(0,3):
                            #curlq = (curlq + M[qi][qj]*MX[qi][qj]*lciv[qi][0][ql])
                            #curlq = (curlq + M[qi][qj]*MZ[qi][qj]*lciv[qi][2][ql])
                            curlq += (M[qi][qj]*MX[ql][qi]*lciv[qj][2][ql])
                            curlq += (M[qi][qj]*MZ[ql][qi]*lciv[qj][0][ql])
                            #curlq += (M[qi][ql]*MX[ql][qi]*lciv[qi][2][ql])
                            #curlq += (M[qi][ql]*MZ[ql][qi]*lciv[qi][0][ql])
                            #print(curlq)
                fcurl.write("%f\n" % (curlq))
        #print(t)

    fcurl.close()
