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
plt.title('Numerical and analytical solution of u(x)',size=18)
plt.plot(x,u)
plt.plot(x,cf)
plt.xlabel('x',size=18)
plt.ylabel('u',size=18)
plt.legend(['Numerical','Analytical'],prop={'size': 18})
plt.show()
