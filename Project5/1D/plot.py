import numpy as np
import matplotlib.pyplot as plt

Euler = False
Back = False
CN = False
Analytical = False
anim1D = True

j = 197
j2 = 1

if anim1D == True:
    u = np.loadtxt("FWEuler.txt")
    uAnaly = np.loadtxt("Analytical1D.txt")

    n = len(u[0,:])
    x = np.linspace(0, 1, n)
    # Put mmatplotlib in interactive mode for animation
    plt.ion()

    # Setup the figure before starting animation
    fig = plt.figure() # Create window
    ax = fig.add_subplot(111) # Add axes
    line, = ax.plot(x, u[0,:], label='u(x)' ) # Fetch the line object
    line2, = ax.plot(x, uAnaly[0,:], label='Analytical')
    ax.set_ylim([0,1])

    for i in range(n):
    #    print len(t), len(u[:,i])
    #    ax.plot(x, u[i,:])
        line.set_ydata(u[i,:]) # Update the y values
        line2.set_ydata(uAnaly[i,:])
        plt.draw() # Update the plot
        plt.pause(0.01)

    # Turn off interactive mode
    plt.ioff()

    # Add show so that windows do not automatically close
    plt.show()

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
