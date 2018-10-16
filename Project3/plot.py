import numpy as np
import matplotlib.pyplot as plt
import sys

infile = sys.argv[1]

x = []
y = []

for lines in open(infile, "r"):
    col = lines.split()
    x.append(float(col[0]))
    y.append(float(col[1]))

plt.plot(x,y)
plt.show()
