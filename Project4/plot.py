import matplotlib.pyplot as plt
import numpy as np
import os

accepts = False
acceptOfT = False
EM = True

file = "data/L20-T1-mcc1e5.txt"

# C_V = np.loadtxt(file, usecols=2, skiprows=2)
# chi = np.loadtxt(file, usecols=3, skiprows=2)
# Mabs = np.loadtxt(file, usecols=4, skiprows=2)
# temp = np.loadtxt("data/test.txt", usecols=4)

if EM == True:
        E = np.loadtxt(file, usecols=0)
        M = np.loadtxt(file, usecols=1)
        mcs = np.linspace(0,len(M),len(M))

        plt.figure()
        plt.plot(mcs,E)
        plt.legend(["Energy"], prop={'size': 15})
        plt.title("Energy with T = 1 for 20x20 lattice", size=15)
        plt.xlabel("Monte Carlo cycles", size=15); plt.ylabel("E", size=15)
        plt.show()

        plt.figure()
        plt.plot(mcs,M)
        plt.legend(["Magnetisation"], prop={'size': 15})
        plt.title("Magnetisation with T = 1 for 20x20 lattice", size=15)
        plt.xlabel("Monte Carlo cycles", size=15); plt.ylabel("M", size=15)
        plt.show()

if accepts == True:
    accepts = np.loadtxt(file, usecols=5)
    mcs = np.linspace(0,len(accepts),len(accepts))
    plt.plot(mcs,accepts)
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



#os.remove("data/test.txt")
