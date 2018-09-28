import matplotlib.pyplot as plt
import sys

infile = sys.argv[1]

eigvals = []

for lines in open(infile, "r"):
    col = lines.split()
    eigvals.append(float(col[0]))


plt.plot(eigvals)
plt.show()
