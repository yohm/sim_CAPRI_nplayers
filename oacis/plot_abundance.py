import pprint
import glob
import re
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

def make_summarized_plot(datfile, names):
    plt.clf()
    dat = np.loadtxt(datfile)
    plt.xlabel('benefit')
    plt.ylabel('abundance')
    plt.ylim(0.0, 1.0)
    legends = ["m=1 defensible", "m=1 efficient", "m=1 other"]
    plot_data = [np.zeros_like(dat[:,1]), np.zeros_like(dat[:,1]), np.zeros_like(dat[:,1])]
    for i,name in enumerate(names):
        if re.search(r'/\d_D$', name):
            plot_data[0] += dat[:,i+1]
        elif re.search(r'/\d_E$', name):
            plot_data[1] += dat[:,i+1]
        elif re.search(r'/\d$', name):
            plot_data[2] += dat[:,i+1]
        else:
            plot_data.append(dat[:,i+1])
            legends.append(name)
    for x,n in zip(plot_data, legends):
        plt.plot(dat[:,0], x, label=n)
    plt.legend()
    plt.savefig(datfile.replace(".dat", ".S.png"))


for d in dat_files:
    make_plot(d, names)
    make_summarized_plot(d, names)

