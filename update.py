import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize

def initupdate(i,Lz,Lt,ks,kt,a,b,c,alpha,beta,gamma,Q3,Q5,z,t,alphatemp,betatemp,gammatemp,Q3temp,Q5temp,coza,cozb):
    sigma = 0.0
    fq3[z,t] = Decimal(-2.0*b*Q5[z,t]*gamma[z,t] + 2.0*Q3[z,t]*(-a + b*alpha[z,t]\
                - b*beta[z,t] + 2.0*c*(Q3[z,t]**2 + Q5[z,t]**2 + 3*alpha[z,t]**2 + beta[z,t]**2 + gamma[z,t]**2)))
    fq5[z,t] = Decimal(-2.0*b*Q3[z,t]*gamma[z,t] + 2.0*Q5[z,t]*(-a + b*(alpha[z,t] + beta[z,t])\
                + 2.0*c*(Q3[z,t]**2 + Q5[z,t]**2 + 3*alpha[z,t]**2 + beta[z,t]**2 + gamma[z,t]**2)))

    falpha[z,t] = Decimal((Q3[z,t]**2 + Q5[z,t]**2)*(b + 12.0*c*alpha[z,t])\
                - 2.0*(-3.0*b*alpha[z,t]**2 - 18.0*c*alpha[z,t]**3 + b*(beta[z,t]**2\
                + gamma[z,t]**2) + 3.0*alpha[z,t]*(a - 2.0*c*(beta[z,t]**2 + gamma[z,t]**2))))

    fbeta[z,t] = Decimal(-(Q3[z,t]**2)*(b - 4.0*c*beta[z,t]) + Q5[z,t]**2*(b + 4.0*c*beta[z,t])\
                - 2.0*beta[z,t]*(a + 2.0*b*alpha[z,t]\
                - 6*c*alpha[z,t]**2 - 2*c*(beta[z,t]**2 + gamma[z,t]**2)))

    fgamma[z,t] = Decimal(-(b*Q3[z,t]*Q5[z,t]) + 2*c*Q3[z,t]**2*gamma[z,t] - gamma[z,t]*(a + 2*b*alpha[z,t]\
                    - 2*c*(Q5[z,t]**2 + 3*alpha[z,t]**2 + beta[z,t]**2 + gamma[z,t]**2)))
    #print("after\n",i,z,falpha[z,t],fbeta[z,t],fgamma[z,t],fq3[z,t],fq5[z,t])

    Q5temp[z,t] = Decimal((-2.0*dz*fq5[z,t] + 4*kt*q0*Q3[cozb,t] - 4*kt*q0*Q3[coza,t] + dz*ks*Q5[cozb,t] + dz*kt*Q5[cozb,t] + \
                dz*ks*Q5[coza,t] + dz*kt*Q5[coza,t])/(2.0*dz*(ks + kt + 2.0*sigma)))

    Q3temp[z,t] =  Decimal((-2.0*dz*fq3[z,t] + dz*ks*Q3[cozb,t] + dz*kt*Q3[cozb,t] + dz*ks*Q3[coza,t] + dz*kt*Q3[coza,t] - \
                4*kt*q0*Q5[cozb,t] + 4*kt*q0*Q5[coza,t])/(2.0*dz*(ks + kt + 2*sigma)))

    gammatemp[z,t] = Decimal((-2.0*dz*fgamma[z,t] + 4*kt*q0*beta[cozb,t] - 4*kt*q0*beta[coza,t] + kt*gamma[cozb,t]\
                    + kt*gamma[coza,t])/(2.0*(kt + dz*sigma)))

    betatemp[z,t] = Decimal((-(dz*fbeta[z,t]) + kt*beta[cozb,t] + kt*beta[coza,t]\
                    - 4.0*kt*q0*gamma[cozb,t] + 4*kt*q0*gamma[coza,t])/(2.0*(kt + dz*sigma)))

    alphatemp[z,t] = Decimal((-(dz*falpha[z,t]) + (2*ks + kt)*(alpha[cozb,t] + alpha[coza,t]))/(2.0*(2*ks + kt + 3*dz*sigma)))

    return(alpha,beta,gamma,Q3,Q5,alphatemp,betatemp,gammatemp,Q3temp,Q5temp)
