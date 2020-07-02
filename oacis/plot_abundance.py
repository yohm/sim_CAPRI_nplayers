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

#pprint.pprint(names)

types = ["m=1 defensible", "m=1 efficient", "m=1 other"]
def species_type(name):
    if re.search(r'/\d_D$', name):
        return 0
    elif re.search(r'/\d_E$', name):
        return 1
    elif re.search(r'/\d$', name):
        return 2
    else:
        return -1

for n in names:
    if species_type(n) < 0:
        types.append(n)

with open("species_types.txt", mode='w') as f:
    f.write('\n'.join(types))

pprint.pprint(types)

def _get_N(fname):
    return int(re.search(r'\d+', fname).group(0))

dat_files = sorted(glob.glob("abundance_*.dat"), key = _get_N)
#pprint.pprint(dat_files)

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

def make_summarized_plot(datfile, names, types):
    plt.clf()
    dat = np.loadtxt(datfile)
    plt.xlabel('benefit')
    plt.ylabel('abundance')
    plt.ylim(0.0, 1.0)
    betas = dat[:,0]
    plot_data = [np.zeros_like(dat[:,1]) for t in types]
    for i,name in enumerate(names):
        t = species_type(name)
        if t >= 0:
            plot_data[t] += dat[:, i+1]
        else:
            idx = types.index(name)
            plot_data[idx] = dat[:,i+1]
    for x,t in zip(plot_data, types):
        plt.plot(betas, x, label=t)
    plt.legend()
    plt.savefig('S_' + datfile.replace(".dat", ".png"))
    d = np.array([betas] + plot_data).T
    return d


s_data = []
ns = [_get_N(fname) for fname in dat_files]

for d in dat_files:
    make_plot(d, names)
    sd = make_summarized_plot(d, names, types)
    np.savetxt('S_' + d, sd)
    s_data.append(sd)

def make_N_dependency_plot(beta_idx, s_data, types, fname):
    plt.clf()
    n_dep = s_data[:,beta_idx,:]
    n_dep[:,0] = ns
    pprint.pprint(n_dep)
    for i,t in enumerate(types):
        plt.plot(ns, n_dep[:,i+1], label=t)
    plt.legend()
    plt.savefig(fname+".png")
    np.savetxt(fname+".dat", n_dep)

s_data = np.array(s_data)
for idx,fname in zip([19,39,59],['N_dep_beta_2','N_dep_beta_3','N_dep_beta_4']):
    make_N_dependency_plot(idx, s_data, types, fname)
