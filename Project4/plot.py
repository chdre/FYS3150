import matplotlib.pyplot as plt
import numpy as np
import os
import glob

accepts = False
acceptOfT = False
EM = False
PE = False
vary = True

file = "data/L100-1e5.txt"
#file2 = "data/T2_4-L20-mcc1e6.txt"


if EM == True:
    E = np.loadtxt(file, usecols=0)
    M = np.loadtxt(file, usecols=1)
    E2 = np.loadtxt(file2, usecols=0)
    M2 = np.loadtxt(file2, usecols=1)
    mcc = np.linspace(0,len(M),len(M))
    mcc2 = np.linspace(0,len(M2), len(M2))

    plt.figure()
    #plt.plot(mcc,E)
    plt.plot(mcc,E2)
    plt.legend(["T=2.4"], prop={'size': 15})
    plt.title("Energy of random state for a 20x20 lattice", size=15)
    plt.xlabel("Monte Carlo cycles", size=15); plt.ylabel("$\\langle E \\rangle/L^2$", size=15)
    plt.show()

    plt.figure()
    #plt.plot(mcc,M)
    plt.plot(mcc,M2)
    plt.legend(["T=2.4"], prop={'size': 15})
    plt.title("Magnetisation of random state for a 20x20 lattice", size=15)
    plt.xlabel("Monte Carlo cycles", size=15); plt.ylabel("$\\langle |M| \\rangle/L^2$", size=15)
    plt.show()

if PE == True:
    E_sorted = np.sort(np.loadtxt(file2, usecols=0))
    mcc = np.linspace(0,len(E_sorted),len(E_sorted))
    bins = 0
    for i in range(0,len(E_sorted)-1):
            if E_sorted[i+1] != E_sorted[i]:
                bins += 1
    plt.hist(E_sorted,bins+10)
    plt.legend(["T=2.4"], prop={'size':15})
    plt.xlabel("Energy", size=15); plt.ylabel("Probability",size=15)
    plt.title("Probability distribution of energy for 20x20 lattice", size=15)
    plt.show()


if accepts == True:
    accepts = np.loadtxt(file, usecols=5)
    accepts2 = np.loadtxt(file2, usecols=5)

    mcc = np.linspace(0,len(accepts),len(accepts))
    plt.loglog(mcc,accepts)
    plt.loglog(mcc,accepts2)
    plt.legend(["T=1", "T=2"], prop={'size': 15})
    plt.title("Accepts in Monte Carlo cycle for 20x20 lattice", size=15)
    plt.xlabel("log(Monte Carlo cycles)", size=15); plt.ylabel("log(Accepts)", size=15)
    plt.show()

if acceptOfT == True:
    acceptT1 = np.loadtxt("data/T1-L20-mcc1e3.txt", usecols=5)
    acceptT2_4 = np.loadtxt("data/T2_4-L20-mcc1e3.txt", usecols=5)
    acceptT12 = np.loadtxt("data/T1-L20-mcc1e4.txt", usecols=5)
    acceptT2_42 = np.loadtxt("data/T2_4-L20-mcc1e4.txt", usecols=5)
    acceptT13 = np.loadtxt("data/T1-L20-mcc1e5.txt", usecols=5)
    acceptT2_43 = np.loadtxt("data/T2_4-L20-mcc1e5.txt", usecols=5)
    acceptT14 = np.loadtxt("data/T1-L20-mcc1e6.txt", usecols=5)
    acceptT2_44 = np.loadtxt("data/T2_4-L20-mcc1e6.txt", usecols=5)

    plt.loglog(1,acceptT1[-1], "ro")
    plt.loglog(1,acceptT12[-1], "bo")
    plt.loglog(1,acceptT13[-1], "co")
    plt.loglog(1,acceptT14[-1], "C4o")
    plt.loglog(2.4, acceptT2_4[-1], "ro")
    plt.loglog(2.4, acceptT2_42[-1], "bo")
    plt.loglog(2.4, acceptT2_43[-1], "co")
    plt.loglog(2.4, acceptT2_44[-1], "C4o")
    plt.legend(["cycles=1e3", "cycles=1e4", "cycles=1e5", "cycles=1e5"], prop={'size': 15})
    plt.title("Accepts states for a 20x20 lattice", size=15)
    plt.xlabel("log(T)", size=15); plt.ylabel("log(Accepts)", size=15)
    plt.show()


if vary == True:
    #files = glob.glob("TempVary-L40-*.txt")
    #for file in files:
    E = np.loadtxt(file, usecols=0)
    Mabs = np.loadtxt(file, usecols=1)
    C_V = np.loadtxt(file, usecols=2)
    chi = np.loadtxt(file, usecols=3)
    T = np.loadtxt(file, usecols=4)

    #plt.figure()
    #for energy in
    plt.plot(T,E)
    plt.legend(["Energy", "M", "C_V", "chi"])
    plt.show()

    plt.plot(T,Mabs)
    plt.show()

    plt.plot(T,C_V)
    plt.show()

    plt.plot(T,chi)
    plt.show()














#os.remove("data/test.txt")
