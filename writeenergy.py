def writeenergy(guess,sig,Lz,Lt,ks,kt,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,s):
    #print("Energy Called")

    dz = 1
    dt = 1
    Qt1,Qt2,Qt3,Qt4,Qt5= (2.0*s/3.0),0.0,0.0,-s/3.0,0.0
    # splay,twist,bend,surface,bulk = 0.0,0.0,0.0,0.0,0.0
    # timeen = 0.0
    enarrayfinal = []

    for t in range(0,Lt):
        for z in range(0,Lz):
            Q1[z,t] = guess[5*t*Lz + 5*z]
            Q2[z,t] = guess[5*t*Lz + 5*z + 1]
            Q3[z,t] = guess[5*t*Lz + 5*z + 2]
            Q4[z,t] = guess[5*t*Lz + 5*z + 3]
            Q5[z,t] = guess[5*t*Lz + 5*z + 4]
    enarrayinit = []

    for t in range(0,Lt):
        splayen,twisten,benden,surfaceen,bulken, timeenen = 0.0,0.0,0.0,0.0,0.0,0.0
        for z in range(0,Lz):
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

        enarrayinit.append(bulken + splayen + twisten + timeenen + surfaceen)

    #enarrayinit = np.array(enarrayinit)
    print(np.shape(enarrayinit))

    energy = (bulk + splay + twist + timeen + surface)# / (Lz*Lt)

    calculategrad.calcgrad(guess,original,GradE,sig,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk)
    np.savetxt("results/ks%ekt%eabc%e%e%esig%d/energyevolution.dat" % (ks,kt,a,b,c,sig), enarrayinit)

    return(energy,GradE)
