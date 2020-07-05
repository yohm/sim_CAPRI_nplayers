import pprint
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


dat = np.loadtxt("cooperation_level.dat")
n_cols = dat.shape[1] - 1
for i in range(n_cols):
    n = i + 2  # N: [3,...]
    col = i + 1
    plt.plot(dat[:,0], dat[:,col], label=n)
plt.legend()
plt.xlabel("benefit")
plt.ylabel("cooperation level")
plt.ylim(0, 1)
plt.savefig("cooperation_level.png")

