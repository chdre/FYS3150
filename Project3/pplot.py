import matplotlib.pyplot as plt
import numpy as np
import sys

data = sys.argv[1]

if data == "momentum":
    angularMom = np.loadtxt("data/test.txt", usecols=0)
    time = np.loadtxt("data/test.txt", usecols=1)


    plt.figure()
    plt.plot(time,angularMom)
    plt.show()

if data == "pos3d":
    x = np.loadtxt("data/3d-escape.txt", usecols=0)
    y = np.loadtxt("data/3d-escape.txt", usecols=1)

    plt.figure()
    plt.plot(x,y)
    plt.title("Earth orbit around sun using Velocity Verlet",size=15)
    plt.xlabel("x",size=15); plt.ylabel("y",size=15)
    plt.legend(['Initial $v = 1.41 * 2\\pi$'],prop={'size': 15})
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


if data == "pos":
    if sys.argv[2] == "verlet":
        x = np.loadtxt("data/3c-verlet.txt", usecols=0)
        y = np.loadtxt("data/3c-verlet.txt", usecols=1)

        plt.figure()
        plt.plot(x,y)
        plt.title("Earth orbit around sun using Velocity Verlet",size=15)
        plt.xlabel("x",size=15); plt.ylabel("y",size=15)
        plt.legend(['Eearth'],prop={'size': 15})
        plt.show()
    if sys.argv[2] == "euler-cromer":
        x = np.loadtxt("data/3c-eulercromer.txt", usecols=0)
        y = np.loadtxt("data/3c-eulercromer.txt", usecols=1)

        plt.figure()
        plt.plot(x,y)
        plt.title("Earth orbit around sun using Euler-Cromer",size=15)
        plt.xlabel("x",size=15); plt.ylabel("y",size=15)
        plt.legend(['Eearth'],prop={'size': 15})
        plt.show()

    if sys.argv[2] == "euler":
        x = np.loadtxt("data/3c-euler.txt", usecols=0)
        y = np.loadtxt("data/3c-euler.txt", usecols=1)

        plt.figure()
        plt.plot(x,y)
        plt.title("Earth orbit around sun using Euler",size=15)
        plt.xlabel("x",size=15); plt.ylabel("y",size=15)
        plt.legend(['Eearth'],prop={'size': 15})
        plt.show()
