import matplotlib.pyplot as plt
import sys

infile = sys.argv[1]

x = []
u = []

for lines in open(infile, "r"):
    col = lines.split()
    x.append(float(col[0]))
    u.append(float(col[1]))

plt.figure()
plt.plot(x,u)
plt.title("LU decomposed u")
plt.xlabel('$x$');plt.ylabel('$u$')
plt.legend(['$u$'])
plt.show()
