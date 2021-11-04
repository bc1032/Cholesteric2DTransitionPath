import math
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import os
import analyse
import energy as en
from matplotlib import colors

os.chdir('/home/bc1032/Desktop/Work/Cholesterics/2DLBFGS/BlueCrystal')


lx,lz,lt = 150,150,200

sig = 1e5
# #Elastic constants
ks = 1e2
kt = 1e2
# #Bulk constants
a=3.0e0
b=2.0e0
c=1.0e0
enbiarr,enunarr, ksarr = [],[],[]
count = -14

def calcen(sig,ks,kt,a,b,c):
    sig = 1e5
    # #Elastic constants
    ks = 1e2
    kt = 1e2
    # #Bulk constants
    a=3.0e0
    b=2.0e0
    c=1.0e0
    guess = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyarray.dat" % (ks,kt,a,b,c,sig))
    analyse.analysis(lx,lz,lt,ks,kt,a,b,c,sig)
    #analyse.writeenergy(guess,sig,lz,lt,ks,kt,a,b,c)
    #while ks <= 1e8:
    #en = np.loadtxt("results/Uniaxial/ks%ekt%eabc%e%e%esig%d/energyevolution.dat" % (ks,kt,a,b,c,sig))
        #kt = ks
        #print(ks)

        # guess = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyarray.dat" % (ks,kt,a,b,c,sig))
        # analyse.analysis(lz,lt,ks,kt,a,b,c,sig)
        # analyse.writeenergy(guess,sig,lz,lt,ks,kt,a,b,c)
        #ks *= 10
        #kt *= 10

def plotaskskt(sig,ks,kt,a,b,c):
    #
    # sig = 1e5
    # #Elastic constants
    # ks = 1e-2
    # kt = 1e-2
    # #Bulk constants
    # a=3.0e0
    # b=2.0e0
    # c=1.0e0

    c1='#1f77b4' #blue
    c2='red' #green
    n=30
    colors = plt.cm.RdYlBu(np.linspace(0,1,n))
    count = -14
    i,x = -14, 0
    enbiarr,power=[],[]
    while i <= 6:

        print(i)
        ks = 1*10**i
        kt = 1*10**i

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/twistevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[1:lt-1]))
        power.append(i)
        plt.plot(enbi, label = '$10^{%d}$' % i, color=colors[i])
        #plt.plot(enun)
        plt.legend(loc='upper left',ncol=2,fontsize='x-small')
        #plt.legend(, ncol=2,handleheight=2.4, labelspacing=0.05, bbox_to_anchor=(1, 0.5))
        #plt.legend()
        plt.xlabel('$t$')

        plt.ylabel('$E_T$')

        i+=1
        x+=1

    plt.savefig("twistevolution.pdf")
    plt.close()
    plt.plot(power,enbiarr)
    #plt.yscale("symlog")
    plt.savefig("maxtwist.pdf")
    plt.close()
    #

    count = -14
    i,x = -14, 0
    enbiarr,power=[],[]
    while i <= 6:

        print(i)
        ks = 1*10**i
        kt = 1*10**i

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/splayevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[2:lt-2]))
        power.append(i)
        plt.plot(enbi, label = '$10^{%d}$' % i, color=colors[i])
        #plt.plot(enun)
        plt.legend(loc='upper left',ncol=2,fontsize='x-small')
        plt.xlabel('$t$')

        plt.ylabel('$E$')

        i+=1
    plt.savefig("splayevolution.pdf")
    plt.close()
    plt.plot(power,enbiarr)
    #plt.yscale("symlog")
    plt.savefig("maxsplaysmall.pdf")
    plt.close()

    count = -14
    i,x = -14, 0
    enbiarr,power=[],[]
    while i <= 6:

        print(i)
        ks = 1*10**i
        kt = 1*10**i

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[2:lt-2]))
        power.append(i)
        plt.plot(enbi, label = '$10^{%d}$' % i, color=colors[i])
        #plt.plot(enun)
        plt.legend(loc='upper left',ncol=2,fontsize='x-small')
        plt.xlabel('$t$')
        plt.yscale("symlog")

        plt.ylabel('$E$')

        i+=1
    plt.savefig("energyevolutionlarge.pdf")
    plt.close()
    plt.plot(power,enbiarr)
    #plt.yscale("symlog")
    plt.savefig("energylarge.pdf")
    plt.close()

    count = -14
    i,x = -14, 0
    enbiarr,power=[],[]
    while i <= 6:

        print(i)
        ks = 1*10**i
        kt = 1*10**i

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/bulkevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[1:lt-1]))
        power.append(i)
        plt.plot(enbi, label = '$10^{%d}$' % i, color=colors[i])
        #plt.plot(enun)
        plt.legend(loc='upper left',ncol=2,fontsize='x-small')
        #plt.legend(, ncol=2,handleheight=2.4, labelspacing=0.05, bbox_to_anchor=(1, 0.5))
        #plt.legend()
        plt.xlabel('$t$')

        plt.ylabel('$E_T$')

        i+=1
        x+=1
    plt.savefig("bulkevolution.pdf")
    plt.close()
    plt.plot(power,enbiarr)
    plt.yscale("symlog")
    plt.savefig("maxbulk.pdf")
    plt.close()

    count = -14
    i,x = -14, 0
    enbiarr,power=[],[]
    while i <= 6:

        print(i)
        ks = 1*10**i
        kt = 1*10**i

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[10:lt-10]))
        plt.plot(enbi)
        #plt.plot(enun)
        plt.xlabel('$t$')

        plt.ylabel('$E$')
        plt.savefig("ks%ekt%eabc%e%e%esig%d/energyevolution.pdf" % (ks,kt,a,b,c,sig))
        plt.close()

        twist = np.loadtxt("ks%ekt%eabc%e%e%esig%d/twistevolution.dat" % (ks,kt,a,b,c,sig))
        plt.plot(twist, label = 'i = %d' % i )
        plt.xlabel('$t$')
        plt.ylabel('$E_{Twist}$')
        plt.legend(loc="upper left")

        plt.savefig("ks%ekt%eabc%e%e%esig%d/twistenergy.pdf" % (ks,kt,a,b,c,sig))
        plt.close()

        splay = np.loadtxt("ks%ekt%eabc%e%e%esig%d/splayevolution.dat" % (ks,kt,a,b,c,sig))
        plt.plot(splay, label = 'i = %d' % i )
        plt.xlabel('$t$')
        plt.ylabel('$E_{Splay}$')
        plt.legend(loc="upper left")

        plt.savefig("ks%ekt%eabc%e%e%esig%d/splayenergy.pdf" % (ks,kt,a,b,c,sig))
        plt.close()

        bulk = np.loadtxt("ks%ekt%eabc%e%e%esig%d/bulkevolution.dat" % (ks,kt,a,b,c,sig))
        plt.plot(bulk, label = 'i = %d' % i )
        plt.xlabel('$t$')
        plt.legend(loc="upper left")

        plt.savefig("ks%ekt%eabc%e%e%esig%d/bulkenergy.pdf" % (ks,kt,a,b,c,sig))
        plt.close()

        i+=1

    return(0)

def calcenvarybulk(sig,ks,kt,a,b,c):
    i = 4
    while a <= 3e12:
        print(i)
        i+=1
        en = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyevolution.dat" % (ks,kt,a,b,c,sig))
        guess = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyarray.dat" % (ks,kt,a,b,c,sig))
        analyse.analysis(lz,lt,ks,kt,a,b,c,sig)
        analyse.writeenergy(guess,sig,lz,lt,ks,kt,a,b,c)
        a*=10
        b*=10
        c*=10
    return(0)

def plotasabc(sig,ks,kt,a,b,c):

    c1='#1f77b4' #blue
    c2='red' #green
    n=30
    colors = plt.cm.RdYlBu(np.linspace(0,1,n))

    count = -6
    i,x = -6, 0
    enbiarr,power=[],[]
    while i <= 12:

        print(i)
        a=3*10**i
        b=2*10**i
        c=1*10**i

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/twistevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[1:lt-1]))
        power.append(i)
        plt.plot(enbi, label = '$10^{%d}$' % i, color=colors[i])
        #plt.plot(enun)
        plt.legend(loc='upper left',ncol=2,fontsize='x-small')
        #plt.legend(, ncol=2,handleheight=2.4, labelspacing=0.05, bbox_to_anchor=(1, 0.5))
        #plt.legend()
        plt.xlabel('$t$')

        plt.ylabel('$E_T$')

        i+=1
        x+=1

    plt.savefig("twistevolution.pdf")
    plt.close()
    plt.plot(power,enbiarr)
    #plt.yscale("symlog")
    plt.savefig("maxtwist.pdf")
    plt.close()

    i = -6
    enbiarr,power=[],[]
    n=20
    colors = plt.cm.RdYlBu(np.linspace(0,1,n))
    while i <= 12:
        print(i)
        a=3*10**i
        b=2*10**i
        c=1*10**i

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/splayevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[2:lt-2]))
        power.append(i)
        plt.plot(enbi, label = '$10^{%d}$' % i, color=colors[i])
        #plt.plot(enun)
        plt.legend(loc='upper left',ncol=2,fontsize='x-small')
        plt.xlabel('$t$')

        plt.ylabel('$E$')

        i+=1
    plt.savefig("splayevolutionsmall.pdf")
    plt.close()
    plt.plot(power,enbiarr)
    #plt.yscale("symlog")
    plt.savefig("maxsplaysmall.pdf")
    plt.close()

    i = -6
    enbiarr,power=[],[]
    n=20
    colors = plt.cm.RdYlBu(np.linspace(0,1,n))
    while i <= 12:
        print(i)
        a=3*10**i
        b=2*10**i
        c=1*10**i

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[2:lt-2]))
        power.append(i)
        plt.plot(enbi, label = '$10^{%d}$' % i, color=colors[i])
        #plt.plot(enun)
        plt.legend(loc='upper left',ncol=2,fontsize='x-small')
        plt.xlabel('$t$')
        plt.yscale("symlog")

        plt.ylabel('$E$')

        i+=1
    plt.savefig("energyevolutionlarge.pdf")
    plt.close()
    plt.plot(power,enbiarr)
    #plt.yscale("symlog")
    plt.savefig("energylarge.pdf")
    plt.close()


    i,x = -6, 0
    enbiarr,power=[],[]
    while i <= 12:

        print(i)
        a=3*10**i
        b=2*10**i
        c=1*10**i

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/bulkevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[1:lt-1]))
        power.append(i)
        plt.plot(enbi, label = '$10^{%d}$' % i, color=colors[i])
        #plt.plot(enun)
        plt.legend(loc='upper left',ncol=2,fontsize='x-small')
        #plt.legend(, ncol=2,handleheight=2.4, labelspacing=0.05, bbox_to_anchor=(1, 0.5))
        #plt.legend()
        plt.xlabel('$t$')

        plt.ylabel('$E_T$')

        i+=1
        x+=1
    plt.savefig("bulkevolution.pdf")
    plt.close()
    plt.plot(power,enbiarr)
    plt.yscale("symlog")
    plt.savefig("maxbulk.pdf")
    plt.close()

    i = -6
    while i <= 12:
        print(i)
        a=3*10**i
        b=2*10**i
        c=1*10**i
        ks = 1e-2
        kt = 1e-2

        enbi = np.loadtxt("ks%ekt%eabc%e%e%esig%d/energyevolution.dat" % (ks,kt,a,b,c,sig))
        enbiarr.append(np.max(enbi[10:lt-10]))
        plt.plot(enbi)
        #plt.plot(enun)
        plt.xlabel('$t$')

        plt.ylabel('$E$')
        plt.savefig("ks%ekt%eabc%e%e%esig%d/energyevolution.pdf" % (ks,kt,a,b,c,sig))
        plt.close()

        twist = np.loadtxt("ks%ekt%eabc%e%e%esig%d/twistevolution.dat" % (ks,kt,a,b,c,sig))
        plt.plot(twist, label = 'i = %d' % i )
        plt.xlabel('$t$')
        plt.ylabel('$E_{Twist}$')
        plt.legend(loc="upper left")

        plt.savefig("ks%ekt%eabc%e%e%esig%d/twistenergy.pdf" % (ks,kt,a,b,c,sig))
        plt.close()

        splay = np.loadtxt("ks%ekt%eabc%e%e%esig%d/splayevolution.dat" % (ks,kt,a,b,c,sig))
        plt.plot(splay, label = 'i = %d' % i )
        plt.xlabel('$t$')
        plt.ylabel('$E_{Splay}$')
        plt.legend(loc="upper left")

        plt.savefig("ks%ekt%eabc%e%e%esig%d/splayenergy.pdf" % (ks,kt,a,b,c,sig))
        plt.close()

        bulk = np.loadtxt("ks%ekt%eabc%e%e%esig%d/bulkevolution.dat" % (ks,kt,a,b,c,sig))
        plt.plot(bulk, label = 'i = %d' % i )
        plt.xlabel('$t$')
        plt.legend(loc="upper left")

        plt.savefig("ks%ekt%eabc%e%e%esig%d/bulkenergy.pdf" % (ks,kt,a,b,c,sig))
        plt.close()

        i+=1

    return(0)
calcen(sig,ks,kt,a,b,c)
#plotaskskt(sig,ks,kt,a,b,c)
#calcenvarybulk(sig,ks,kt,a,b,c)
#plotasabc(sig,ks,kt,a,b,c)
