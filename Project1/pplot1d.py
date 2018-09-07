import matplotlib.pyplot as plt
import sys

infile = sys.argv[1]

log10h = []
eps = []

for lines in open(infile, "r"):
    col = lines.split()
    eps.append(float(col[0]))
    log10h.append(float(col[1]))

plt.figure()
plt.plot(log10h,eps)
plt.xlabel('$\\log_{10}(h)$');plt.ylabel('$\\epsilon_i$')
plt.legend(['$\\epsilon_i$'])
plt.show()
