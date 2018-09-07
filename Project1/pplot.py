import matplotlib.pyplot as plt
import sys

infile = sys.argv[1]

x = []
u = []
cf = []

for lines in open(infile, "r"):
    col = lines.split()
    x.append(float(col[0]))
    u.append(float(col[1]))
    cf.append(float(col[2]))

plt.figure()
plt.plot(x,u)
plt.plot(x,cf)
plt.legend(['u','cf'])

plt.show()
