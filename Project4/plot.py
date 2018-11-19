import matplotlib.pyplot as plt
import numpy as np
import os

accepts = False
acceptOfT = False
EM = False
PE = False
vary = False
variance = False

if EM == True:
    file = "data/L20-GS-T1-1e6.txt"
    file2 = "data/L20-GS-T2_4-1e6.txt"

    E = np.loadtxt(file, usecols=0)
    M = np.loadtxt(file, usecols=1)
    E2 = np.loadtxt(file2, usecols=0)
    M2 = np.loadtxt(file2, usecols=1)
    mcc = np.linspace(0,len(M),len(M))

    plt.figure()
    #plt.plot(mcc,E)
    plt.plot(mcc,E2)
    plt.legend(["T=1"], prop={'size': 15})
    plt.title("Energy for a 20x20 lattice (ground state)", size=15)
    plt.xlabel("Monte Carlo cycles", size=15); plt.ylabel("$ E /L^2$", size=15)
    plt.show()

    plt.figure()
    #plt.plot(mcc,M)
    plt.plot(mcc,M2)
    plt.legend(["T=1"], prop={'size': 15})
    plt.title("Magnetisation for a 20x20 lattice (ground state)", size=15)
    plt.xlabel("Monte Carlo cycles", size=15); plt.ylabel("$|M|/L^2$", size=15)
    plt.show()

if PE == True:
    E_sorted = np.sort(np.loadtxt(dfile, usecols=0))
    mcc = np.linspace(0,len(E_sorted),len(E_sorted))
    bins = 0
    for i in range(0,len(E_sorted)-1):
            if E_sorted[i+1] != E_sorted[i]:
                bins += 1
    plt.hist(E_sorted,bins+10)
    plt.legend(["T=1"], prop={'size':15})
    plt.xlabel("Energy", size=15); plt.ylabel("Probability",size=15)
    plt.title("Probability distribution of energy for 20x20 lattice", size=15)
    plt.show()


if accepts == True:
    file = "data/L20-T1-1e5.txt"
    file2 = "data/L20-T2_4-1e5.txt"
    accepts = np.loadtxt(file, usecols=5)
    accepts2 = np.loadtxt(file2, usecols=5)

    mcc = np.linspace(0,len(accepts),len(accepts))
    plt.loglog(mcc,accepts)
    plt.loglog(mcc,accepts2)
    plt.legend(["T=1", "T=2.4"], prop={'size': 15})
    plt.title("Accepts in Monte Carlo cycle for 20x20 lattice", size=15)
    plt.xlabel("log(Monte Carlo cycles)", size=15); plt.ylabel("log(Accepts)", size=15)
    plt.show()

if acceptOfT == True:
    acceptT1 = np.loadtxt("data/old/L20-T1-mcc1e3-accepts.txt", usecols=5)
    acceptT2_4 = np.loadtxt("data/old/L20-T2_4-mcc1e3-accepts.txt", usecols=5)
    acceptT12 = np.loadtxt("data/old/L20-T1-mcc1e4-accepts.txt", usecols=5)
    acceptT2_42 = np.loadtxt("data/old/L20-T2_4-mcc1e4-accepts.txt", usecols=5)
    acceptT13 = np.loadtxt("data/old/L20-T1-mcc1e5-accepts.txt", usecols=5)
    acceptT2_43 = np.loadtxt("data/old/L20-T2_4-mcc1e5-accepts.txt", usecols=5)
    acceptT14 = np.loadtxt("data/old/L20-T1-mcc1e6-accepts.txt", usecols=5)
    acceptT2_44 = np.loadtxt("data/old/L20-T2_4-mcc1e6-accepts.txt", usecols=5)

    plt.loglog(1,acceptT1[-1], "ro")
    plt.loglog(1,acceptT12[-1], "bo")
    plt.loglog(1,acceptT13[-1], "co")
    plt.loglog(1,acceptT14[-1], "C4o")
    plt.loglog(2.4, acceptT2_4[-1], "ro")
    plt.loglog(2.4, acceptT2_42[-1], "bo")
    plt.loglog(2.4, acceptT2_43[-1], "co")
    plt.loglog(2.4, acceptT2_44[-1], "C4o")
    plt.legend(["cycles=1e3", "cycles=1e4", "cycles=1e5", "cycles=1e6"], prop={'size': 15})
    plt.title("Accepts states for a 20x20 lattice", size=15)
    plt.xlabel("log(T)", size=15); plt.ylabel("log(Accepts)", size=15)
    plt.show()


if vary == True:
    file = "data/test/L40-1e6-ext.txt"
    file2 = "data/test/L60-1e6-ext.txt"
    file3 = "data/test/L80-1e6-ext.txt"
    file4 = "data/test/L100-1e6-ext.txt"

    E, Mabs, C_V, chi, T = np.loadtxt(file, usecols=(0,1,2,3,4), unpack=True)
    E2, Mabs2, C_V2, chi2 = np.loadtxt(file2, usecols=(0,1,2,3), unpack=True)
    E3, Mabs3, C_V3, chi3 = np.loadtxt(file3, usecols=(0,1,2,3), unpack=True)
    E4, Mabs4, C_V4, chi4 = np.loadtxt(file4, usecols=(0,1,2,3), unpack=True)

    plt.figure()
    plt.plot(T, E, T, E2, T, E3, T, E4)
    plt.legend(["L=40", "L=60", "L=80", "L=100"], prop={'size':15})
    plt.title('Mean energy for various lattice sizes', size=15)
    plt.xlabel('T [kT/J]', size=15); plt.ylabel('$\\langle E \\rangle/L^2$', size=15)
    plt.show()

    plt.figure()
    plt.plot(T, Mabs, T, Mabs2, T, Mabs3, T, Mabs4)
    plt.legend(["L=40", "L=60", "L=80", "L=100"], prop={'size':15})
    plt.title('Mean value of magnetic moment for various lattice sizes', size=15)
    plt.xlabel('T [kT/J]', size=15); plt.ylabel('$\\langle|M|\\rangle/L^2$', size=15)
    plt.show()

    plt.figure()
    plt.plot(T, C_V, T, C_V2, T, C_V3, T, C_V4)
    plt.legend(["L=40", "L=60", "L=80", "L=100"], prop={'size':15})
    plt.title('Heat capacity for various lattice sizes', size=15)
    plt.xlabel('T [kT/J]', size=15); plt.ylabel('$C_V/L^2$', size=15)
    plt.show()

    plt.figure()
    plt.plot(T, chi, T, chi2, T, chi3, T, chi4)
    plt.legend(["L=40", "L=60", "L=80", "L=100"], prop={'size':15})
    plt.title('Susceptibility for various lattice sizes', size=15)
    plt.xlabel('T [kT/J]', size=15); plt.ylabel('$\\chi/L^2$', size=15)
    plt.show()

if variance == True:
    file = "data/L20-T1-1e6.txt"
    file2 = "data/L20-T2_4-1e6.txt"

    var = np.loadtxt(file, usecols=4)
    var2 = np.loadtxt(file2, usecols=4)

    mcc = np.linspace(0,len(var),len(var))

    print(var[-1],var2[-1])

    plt.figure()
    plt.plot(mcc, var)
    plt.plot(mcc, var2)
    plt.legend(["T=1", "T=2.4"])
    plt.show()














#os.remove("data/test.txt")
