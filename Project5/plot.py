import numpy as np
import matplotlib.pyplot as plt

Euler = True
Back = True
CN = True

if Euler == True:
    u = np.loadtxt("FWEuler.txt", usecols=10)
    u2 = np.loadtxt("FWEuler.txt", usecols=2)

    plt.plot(u)
    plt.plot(u2)
    plt.show()


if Back == True:
    u = np.loadtxt("BWEuler.txt", usecols=10)
    u2 = np.loadtxt("BWEuler.txt", usecols=2)

    plt.plot(u)
    plt.plot(u2)
    plt.show()

if CN == True:
    u = np.loadtxt("CrankNic.txt", usecols=10)
    u2 = np.loadtxt("CrankNic.txt", usecols=2
    )

    plt.plot(u)
    plt.plot(u2)
    plt.show()
