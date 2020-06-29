import pprint
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


names = []
with open("species.txt") as f:
    names = [a.strip() for a in f.readlines() if len(a.strip()) > 0]

pprint.pprint(names)

dat_files = sorted(glob.glob("abundance_*.dat"))
pprint.pprint(dat_files)

def make_plot(datfile, names):
    plt.clf()
    dat = np.loadtxt(datfile)
    th = 0.05
    plt.xlabel('benefit')
    plt.ylabel('abundance')
    plt.ylim(0.0, 1.0)
    for i,name in enumerate(names):
        if np.max(dat[:,i+1]) > th:
            plt.plot(dat[:,0], dat[:,i+1], label=name)
    plt.legend()
    plt.savefig(datfile.replace(".dat", ".png"))

for d in dat_files:
    make_plot(d, names)

