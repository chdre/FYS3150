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
plt.title('Logarithmic plot of error as a function of step size h',size=18)
plt.plot(log10h,eps)
plt.xlabel('$\\log_{10}(h)$', size=18);plt.ylabel('$\\epsilon_i$',size=18)
plt.legend(['$\\log10(\\epsilon_i)$'],prop={'size': 18})
plt.show()
