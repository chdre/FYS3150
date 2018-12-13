import numpy as np
import matplotlib.pyplot as plt

Euler = False
Back = True
CN = False
Analytical = False
anim1D = False

j = 197
j2 = 1

if Analytical == True:
    u = np.loadtxt("Analytical1D.txt");
    plt.plot(u)
    plt.show()

if Euler == True:
    u = np.loadtxt("FWEuler.txt")

    n = len(u[0,:])
    x = np.linspace(0,1,n)

    plt.plot(x,u[j,:])
    plt.plot(x,u[j2,:])
    plt.show()


if Back == True:
    u = np.loadtxt("BWEuler.txt")
    u = u[:,1:]


    n = len(u[0,:])
    x = np.linspace(0,1,n)

    plt.plot(x,u[j,:])
    plt.plot(x,u[j2,:])
    plt.show()

if CN == True:
    u = np.loadtxt("CrankNic.txt")
    u = u[:,1:]

    n = len(u[0,:])
    x = np.linspace(0,1,n)

    plt.plot(x,u[j,:])
    plt.plot(x,u[j2,:])
    plt.show()
