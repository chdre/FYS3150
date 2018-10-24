import matplotlib.pyplot as plt
import numpy as np
import sys

data = sys.argv[1]

params = {'legend.fontsize':'x-large',
        'axes.labelsize': 'x-large',
        'axes.titlesize': 'x-large',
        'xtick.labelsize': 'x-large',
        'ytick.labelsize': 'x-large'}
plt.rcParams.update(params)

params = {'legend.fontsize': 'x-large',
        'axes.labelsize': 'x-large',
        'axes.titlesize': 'x-large',
        'xtick.labelsize': 'x-large',
        'ytick.labelsize': 'x-large'}
plt.rcParams.update(params)

if data == "momentum":
    angularMom = np.loadtxt("data/test.txt", usecols=0)
    time = np.loadtxt("data/test.txt", usecols=1)


    plt.figure()
    plt.title("Momentum of Earth-Sun system using Velocity Verlet")
    plt.plot(time,angularMom)
    plt.xlabel('Time'); plt.ylabel('p')
    plt.show()

if data == "pos3d":
    x = np.loadtxt("data/3d-escape.txt", usecols=0)
    y = np.loadtxt("data/3d-escape.txt", usecols=1)

    plt.figure()
    plt.plot(x,y)
    plt.title("Earth orbiting the Sun",size=15)
    plt.xlabel("x",size=15); plt.ylabel("y",size=15)
    plt.legend(['Initial $v = 1.41 * 2\\pi$'],prop={'size': 15})
    plt.show()

if data == "pos3dr3":
    x = np.loadtxt("data/3d-newforce.txt", usecols=0)
    y = np.loadtxt("data/3d-newforce.txt", usecols=1)

    plt.figure()
    plt.plot(x,y)
    plt.title("Earth orbiting the Sun using Velocity Verlet",size=15)
    plt.xlabel("x",size=15); plt.ylabel("y",size=15)
    plt.legend(['Earth'],prop={'size': 15})
    plt.show()

if data == "pos3e":
    x = np.loadtxt("data/3e.txt", usecols=0)
    y = np.loadtxt("data/3e.txt", usecols=1)
    x2 = np.loadtxt("data/3e.txt", usecols=2)
    y2 = np.loadtxt("data/3e.txt", usecols=3)


    plt.figure()
    plt.plot(x,y)
    plt.plot(x2,y2)
    plt.title("Earth and Jupiter orbiting the Sun using Velocity Verlet",size=15)
    plt.xlabel("x",size=15); plt.ylabel("y",size=15)
    plt.legend(['Earth','Jupiter'],prop={'size': 15})
    plt.show()

if data == "pos3e2":
    x = np.loadtxt("data/3e2.txt", usecols=0)
    y = np.loadtxt("data/3e2.txt", usecols=1)
    x2 = np.loadtxt("data/3e2.txt", usecols=2)
    y2 = np.loadtxt("data/3e2.txt", usecols=3)


    plt.figure()
    plt.plot(x,y)
    plt.plot(x2,y2)
    plt.plot(0,0,"ko")
    plt.title("Earth and Jupiter orbiting the Sun",size=15)
    plt.xlabel("x",size=15); plt.ylabel("y",size=15)
    plt.legend(['Earth','Jupiter', 'Sun'],prop={'size': 15},loc='upper right')
    plt.show()

if data == "pos3f":
    x = np.loadtxt("data/3f.txt", usecols=0)
    y = np.loadtxt("data/3f.txt", usecols=1)
    x2 = np.loadtxt("data/3f.txt", usecols=2)
    y2 = np.loadtxt("data/3f.txt", usecols=3)
    x3 = np.loadtxt("data/3f.txt", usecols=4)
    y3 = np.loadtxt("data/3f.txt", usecols=5)
    x4 = np.loadtxt("data/3f.txt", usecols=6)
    y4 = np.loadtxt("data/3f.txt", usecols=7)


    plt.figure()
    plt.plot(x,y)
    plt.plot(x2,y2)
    plt.plot(x3,y3)
    plt.plot(x4,y4)
    plt.title("Solar system of Earth, Jupiter and Sun",size=15)
    plt.xlabel("x",size=15); plt.ylabel("y",size=15)
    plt.legend(['Earth','Jupiter', 'Sun', 'Center of Mass'],prop={'size': 15},loc='upper right')
    plt.show()

if data == "pos3g":
    x = np.loadtxt("data/3g.txt", usecols=0)
    y = np.loadtxt("data/3g.txt", usecols=1)
    x2 = np.loadtxt("data/3g.txt", usecols=2)
    y2 = np.loadtxt("data/3g.txt", usecols=3)


    plt.figure()
    plt.plot(x,y)
    plt.plot(x2,y2)
    plt.title("Solar system of Mercury and Sun",size=15)
    plt.xlabel("x",size=15); plt.ylabel("y",size=15)
    plt.legend(['Mercury', 'Sun'],prop={'size': 15},loc='upper right')
    plt.show()

if data == "energy":
    kinetic = np.loadtxt("data/energy.txt", usecols=0)
    potential = np.loadtxt("data/energy.txt", usecols=1)
    time = np.loadtxt("data/energy.txt", usecols=2)

    plt.figure()
    plt.plot(time,kinetic)
    plt.plot(time,potential)

    plt.title("Energy of Earth orbiting the Sun",size=15)
    plt.xlabel("time [years]",size=15); plt.ylabel("Energy",size=15)
    plt.legend(['Kinetic', 'Potential'],prop={'size': 15})
    plt.show()


if data == "pos":
    if sys.argv[2] == "verlet":
        x = np.loadtxt("data/3c-verlet.txt", usecols=0)
        y = np.loadtxt("data/3c-verlet.txt", usecols=1)

        plt.figure()
        plt.plot(x,y)
        plt.title("Earth orbiting the Sun using Velocity Verlet",size=15)
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
        plt.title("Earth orbiting the Sun using Euler",size=15)
        plt.xlabel("x",size=15); plt.ylabel("y",size=15)
        plt.legend(['Eearth'],prop={'size': 15})
        plt.show()
