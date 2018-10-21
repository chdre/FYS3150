import matplotlib.pyplot as plt
import sys

infile = sys.argv[1]

x = []
y = []


for lines in open(infile, "r"):
    col = lines.split()
    x.append(float(col[0]))
    y.append(float(col[1]))

plt.figure()
plt.plot(x,y)
plt.title("Earth orbit around sun using Velocity Verlet",size=15)
plt.xlabel("x",size=15); plt.ylabel("y",size=15)
plt.legend(['Eearth'],prop={'size': 15})
plt.show()
