import matplotlib.pyplot as plt
import numpy as np
import sys

data = sys.argv[1]

if data == "momentum":
    angularMom = np.loadtxt("data/test.txt", usecols=0)
    time = np.loadtxt("data/energy.txt", usecols=1)


    plt.figure()
    plt.plot(time,angularMom)
    plt.show()


if data == "energy":
    kinetic = np.loadtxt("data/energy.txt", usecols=0)
    potential = np.loadtxt("data/energy.txt", usecols=1)
    time = np.loadtxt("data/energy.txt", usecols=2)

    plt.figure()
    plt.plot(time,kinetic)
    plt.plot(time,potential)

    plt.title("Energy of Earth in orbit around sun",size=15)
    plt.xlabel("time [years]",size=15); plt.ylabel("Energy",size=15)
    plt.legend(['Kinetic', 'Potential'],prop={'size': 15})
    plt.show()
