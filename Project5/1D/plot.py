import numpy as np
import matplotlib.pyplot as plt

Euler = True
Back = True
CN = True

j = 103
j2 = 30

if Euler == True:
    u = np.loadtxt("FWEuler.txt")
    u = u[:,1:]

    print(u.shape)
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
