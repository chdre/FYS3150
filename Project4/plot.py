import matplotlib.pyplot as plt
import numpy as np
import os
import glob

accepts = False
acceptOfT = False
EM = False
sort = False
vary = True

file = "data/TempVary-L40.txt"


if EM == True:
    E = np.loadtxt(file, usecols=0)
    M = np.loadtxt(file, usecols=1)
    mcc = np.linspace(0,len(M),len(M))

    plt.figure()
    plt.plot(mcc,E)
    plt.legend(["Energy"], prop={'size': 15})
    plt.title("Energy with T = 1 for 20x20 lattice", size=15)
    plt.xlabel("Monte Carlo cycles", size=15); plt.ylabel("E", size=15)
    plt.show()

    plt.figure()
    plt.plot(mcc,M)
    plt.legend(["Magnetisation"], prop={'size': 15})
    plt.title("Magnetisation with T = 1 for 20x20 lattice", size=15)
    plt.xlabel("Monte Carlo cycles", size=15); plt.ylabel("M", size=15)
    plt.show()

if sort == True:
    E_sorted = np.sort(np.loadtxt(file, usecols=0))
    mcc = np.linspace(0,len(E_sorted),len(E_sorted))

    d = ({1,E_sorted[0]})
    j = 0
    




    plt.plot(E_sorted)
    plt.show()


if accepts == True:
    accepts = np.loadtxt(file, usecols=5)
    mcc = np.linspace(0,len(accepts),len(accepts))
    plt.plot(mcc,accepts)
    plt.legend(["Accepts"], prop={'size': 15})
    plt.title("Accepts in Monte Carlo cycle for T=2.4 and 20x20 lattice", size=15)
    plt.xlabel("Monte Carlo cycles", size=15); plt.ylabel("Accepts", size=15)
    plt.show()

if acceptOfT == True:
    acceptT1 = np.loadtxt("data/L20-T1-mcc1e3-accepts.txt", usecols=5)
    acceptT2_4 = np.loadtxt("data/L20-T2_4-mcc1e3-accepts.txt", usecols=5)
    acceptT12 = np.loadtxt("data/L20-T1-mcc1e4-accepts.txt", usecols=5)
    acceptT2_42 = np.loadtxt("data/L20-T2_4-mcc1e4-accepts.txt", usecols=5)
    acceptT13 = np.loadtxt("data/L20-T1-mcc1e5-accepts.txt", usecols=5)
    acceptT2_43 = np.loadtxt("data/L20-T2_4-mcc1e5-accepts.txt", usecols=5)


    plt.plot(1,acceptT1[-1], "bo")
    plt.plot(1,acceptT12[-1], "go")
    plt.plot(1,acceptT13[-1], "co")
    plt.plot(2.4, acceptT2_42[-1], "go")
    plt.plot(2.4, acceptT2_4[-1], "bo")
    plt.plot(2.4, acceptT2_43[-1], "co")
    plt.legend(["cycles=1e3", "cycles=1e4", "cycles=1e5"], prop={'size': 15})
    plt.title("Zoom of accepts in MC cycle for T=1", size=15)
    plt.xlabel("T [kT/J]", size=15); plt.ylabel("Accepts", size=15)
    plt.show()


if vary == True:
    #files = glob.glob("TempVary-L40-*.txt")
    #for file in files:
    E = np.loadtxt(file, usecols=0)
    Mabs = np.loadtxt(file, usecols=1)
    C_V = np.loadtxt(file, usecols=2)
    chi = np.loadtxt(file, usecols=3)
    T = np.loadtxt(file, usecols=4)

    plt.plot(T,E)
    plt.show()














#os.remove("data/test.txt")
