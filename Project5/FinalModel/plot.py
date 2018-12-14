import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

Euler = False
anim = True

j = 10
j2 = 30

if anim == True:
    u = np.loadtxt("BWEuler2D.txt")
    n = len(u[0,:])
    x = np.linspace(0, 1, n)

    t_steps = int(len(u[:,0])/len(u[0,:]))  # get number of time steps from matrix
    mat = np.zeros((t_steps, n, n))
    # Create matrices for individual time steps
    for t in range(t_steps):
        start = n*t
        end = start+n
        mat[t] = u[start:end]

    fig = plt.figure(figsize=(8,8))
    im = plt.imshow(mat[0], cmap=cm.coolwarm, animated=True)
    #plt.axis('equal')

    #legend
    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('$T(x, y)$', rotation=270, size=15)
    plt.xlabel('$x$', size=15); plt.ylabel('$y$', size=15)
    plt.title('Temperature distribution in a %s x %s grid' % (n,n), size=17)

    i = 0
    def updatefig(*args):
        global i
        print (i)
        if (i < t_steps):
            i += 1
        else:
            i=0         # reset animation
        im.set_array(mat[int(i)-1])
        return im,


    save = True
    ani = animation.FuncAnimation(fig, updatefig, interval=1, blit=True, frames=int(t_steps))

    if save == True:
        # Set up formatting for the movie files
        Writer = animation.writers['ffmpeg']      # requires ffmpeg to be installed
        writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        ani.save('2D_anim.mp4', writer=writer)
    else:
        plt.show()


if Euler == True:
    u = np.loadtxt("FWEuler2D.txt")
    n = len(u[0,:])

    vals = []

    t_steps = int(len(u[:,0])/n)

    index = 0 # dummy index to add matrices to vals array

    for t in range(t_steps):
        start = n*t
        end = start+n
        #print(start, end, t, t_steps)
        #print(u[start:end])

        vals.append(u[start:end])

        index += 1

    vals = np.array(vals)

    print(vals[-1][0], vals[-1][2])

    x = np.linspace(0,1,n)

    plt.plot(vals[-1])
    plt.show()
