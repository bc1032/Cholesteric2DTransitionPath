import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import os
import analyse
import energy as en

os.chdir('/home/bc1032/Desktop/Work/Cholesterics/LBFGS/results/DisparateUpdate/results/')


lz,lt = 100,800

sig = 1e5
#Elastic constants
ks = 1e-14
kt = 1e-14
#Bulk constants
a=3e-4
b=2e-4
c=1e-4
enbiarr,enunarr, ksarr = [],[],[]
count = -14

def calcen():
    ks = 1e-14
    kt = 1e-14
    a=3e-4
    b=2e-4
    c=1e-4
    while ks <= 1e10:
    #en = np.loadtxt("results/Uniaxial/ks%ekt%eabc%e%e%esig%d/energyevolution.dat" % (ks,kt,a,b,c,sig))
        #kt = ks
        print(ks)

        guess = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyarray.dat" % (ks,kt,a,b,c,sig))
        analyse.analysis(lz,lt,ks,kt,a,b,c,sig)
        analyse.writeenergy(guess,sig,lz,lt,ks,kt,a,b,c)
        #ks *= 10
        ks *= 10

def plotaskskt():

    ks = 1e-14
    kt = 1e-14
    count = -14
    a=3e-4
    b=2e-4
    c=1e-4

    while kt <= 1e6:
        print("plot %f" % ks)

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[10:lt-10]))
        plt.plot(enbi)
        #plt.plot(enun)
        plt.xlabel('$t$')

        plt.ylabel('$E$')
        plt.savefig("energyevolutionks%ekt%e.pdf" % (ks,kt))
        #plt.show()
        plt.close()
        ksarr.append(count)

        count += 1
        #ks *= 10
        kt *= 10

    plt.plot(ksarr, enbiarr)
    #plt.plot(ksarr, enunarr)
    plt.xlabel('$k_s, k_t$')
    plt.ylabel('$E_{max}$')
    plt.savefig("maxenergyasratio.pdf")
    #plt.show()
    plt.close()

    ks = 1e-14
    kt = 1e-14

    while ks <= 1e6:
        twist = np.loadtxt("ks%ekt%eabc%e%e%esig%d/twistevolution.dat" % (ks,kt,a,b,c,sig))
        plt.plot(twist, label = 'k_T = %e' % kt )
        plt.xlabel('$t$')
        plt.legend(loc="upper left")

        ks *= 6

        #plt.show()
    plt.savefig("twistenergy.pdf")
    plt.close()

    splayarr2,ksarr2=[],[]
    ks = 1e-14
    kt = 1e-14
    while kt <= 1e6:

        splay = np.loadtxt("ks%ekt%eabc%e%e%esig%d/splayevolution.dat" % (ks,kt,a,b,c,sig))
        plt.plot(splay, label = 'k_s = %e' % ks )
        plt.xlabel('$t$')
        plt.legend(loc="upper left")
        splayarr2.append(np.max(splay))
        ksarr2.append(ks)
        kt *= 10

    plt.show()
    plt.savefig("splayenergy.pdf")
    plt.close()
    print(np.shape(splayarr2))
    plt.plot(ksarr2,splayarr2)
    plt.show()
    plt.savefig("maxsplayenergy.pdf")
    plt.close()
    ks = 1e-14
    kt = 1e-14
    while kt <= 1e6:
        bulk = np.loadtxt("ks%ekt%eabc%e%e%esig%d/bulkevolution.dat" % (ks,kt,a,b,c,sig))
        plt.plot(bulk, label = 'k_s = %e' % ks )
        plt.xlabel('$t$')
        plt.legend(loc="upper left")
        #plt.ylabel('$E_{\mathrm{twist}}$')
        #plt.legend('k_T = %e' % kt)
        #ks *= 10
        kt *= 10

        #plt.show()
    plt.savefig("bulkenergy.pdf")
    plt.close()


    return(0)

#calcen()
plotaskskt()
