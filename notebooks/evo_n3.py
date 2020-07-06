# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

plt.rcParams['font.family'] ='sans-serif'
plt.rcParams["figure.subplot.left"] = 0.22
plt.rcParams["figure.subplot.right"] = 0.95
plt.rcParams["figure.subplot.bottom"] = 0.20
plt.rcParams["figure.subplot.top"] = 0.95
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['axes.labelsize'] = 30
plt.rcParams['font.size'] = 22
plt.rcParams['axes.linewidth'] = 1.2

data_dir = "evo_n3_data"

# http://18.177.82.245:3000/runs/5efd76b5e09e252b4d21b91a
dat = np.loadtxt(f"{data_dir}/m1_S_abundance_30.dat")

plt.ylim(-0.02, 1.02)
plt.yticks([0,0.5,1.0])
plt.xlabel("benefit-to-cost ratio")
plt.ylabel("abundance")
plt.plot(dat[:,0], dat[:,1], label="m=1 defensible")
plt.plot(dat[:,0], dat[:,2], label="m=1 efficient")
plt.plot(dat[:,0], dat[:,3], label="m=1 others")
#plt.legend()
plt.savefig(f"{data_dir}/m1_S_abundance_30.pdf")

# num_defensible = 4
# num_efficient = 14
plt.ylim(1.0e-1, 1.0e1)
plt.xlabel("benefit-to-cost ratio")
plt.yscale('log')
#plt.yticks(color="None")
plt.ylabel("Rescaled Abundance")
plt.plot(dat[:,0], dat[:,1]/4 *64.0, label="m=1 defensible")
plt.plot(dat[:,0], dat[:,2]/14*64.0, label="m=1 efficient")
plt.plot(dat[:,0], dat[:,3]/46*64.0, label="m=1 others")

dat = np.loadtxt(f"{data_dir}/m1_N_dep_beta_3.dat")

plt.ylim(-0.02, 1.02)
plt.xlabel("population size")
#plt.xscale('log')
plt.yticks(color="None")
#plt.ylabel("Abundance")
plt.plot(dat[:,0], dat[:,1], label="m=1 defensible")
plt.plot(dat[:,0], dat[:,2], label="m=1 efficient")
plt.plot(dat[:,0], dat[:,3], label="m=1 others")
plt.savefig(f"{data_dir}/m1_N_dep_beta_3.pdf")

# num_defensible = 4
# num_efficient = 14
plt.ylim(1.0e-1, 1.0e1)
plt.xlabel("population size")
plt.yscale('log')
#plt.yticks(color="None")
plt.ylabel("Rescaled Abundance")
plt.plot(dat[:,0], dat[:,1]/4 *64.0, label="m=1 defensible")
plt.plot(dat[:,0], dat[:,2]/14*64.0, label="m=1 efficient")
plt.plot(dat[:,0], dat[:,3]/46*64.0, label="m=1 others")
#plt.savefig("evo/m1_N_dep_beta_3.pdf")


dat = np.loadtxt(f"{data_dir}/m1_coop_level.dat")

plt.ylim(-0.02, 1.02)
plt.xlabel("benefit-to-cost ratio")
plt.ylabel("cooperation level")
plt.plot(dat[:,0], dat[:,28])
plt.savefig(f"{data_dir}/m1_coop_level.pdf")

# http://18.177.82.245:3000/runs/5efd76b5e09e252b4d21b91c
dat = np.loadtxt(f"{data_dir}/CAPRI3_S_abundance_30.dat")

plt.ylim(-0.02, 1.02)
plt.yticks([0,0.5,1.0])
plt.xlabel("benefit-to-cost ratio")
plt.ylabel("abundance")
plt.plot(dat[:,0], dat[:,1], label="m=1 defensible")
plt.plot(dat[:,0], dat[:,2], label="m=1 efficient")
plt.plot(dat[:,0], dat[:,3], label="m=1 others")
plt.plot(dat[:,0], dat[:,4], label="CAPRI-3")
plt.savefig(f"{data_dir}/CAPRI3_S_abundance_30.pdf")

# num_defensible = 4
# num_efficient = 14
plt.ylim(1.0e-2, 1.0e2)
plt.xlabel("benefit-to-cost ratio")
plt.yscale('log')
#plt.yticks(color="None")
plt.ylabel("Rescaled Abundance")
plt.plot(dat[:,0], dat[:,1]/4 *64.0, label="m=1 defensible")
plt.plot(dat[:,0], dat[:,2]/14*64.0, label="m=1 efficient")
plt.plot(dat[:,0], dat[:,3]/46*64.0, label="m=1 others")
plt.plot(dat[:,0], dat[:,4]   *64.0, label="CAPRI-3")

dat = np.loadtxt(f"{data_dir}/CAPRI3_N_dep_beta_3.dat")

plt.ylim(-0.02, 1.02)
plt.xlabel("N")
#plt.xscale('log')
plt.yticks(color="None")
#plt.ylabel("Abundance")
plt.plot(dat[:,0], dat[:,1], label="m=1 defensible")
plt.plot(dat[:,0], dat[:,2], label="m=1 efficient")
plt.plot(dat[:,0], dat[:,3], label="m=1 others")
plt.plot(dat[:,0], dat[:,4], label="CAPRI-3")
plt.savefig(f"{data_dir}/CAPRI3_N_dep_beta_3.pdf")

dat = np.loadtxt(f"{data_dir}/CAPRI3_coop_level.dat")

plt.ylim(-0.02, 1.02)
plt.xlabel("benefit-to-cost ratio")
plt.ylabel("cooperation level")
plt.plot(dat[:,0], dat[:,28])
plt.savefig(f"{data_dir}/CAPRI3_coop_level.pdf")

# http://18.177.82.245:3000/runs/5efda208e09e25f4fb41dc5a
dat = np.loadtxt(f"{data_dir}/ALL_S_abundance_30.dat")

plt.ylim(-0.02, 1.02)
plt.yticks([0,0.5,1.0])
plt.xlabel("benefit-to-cost ratio")
plt.ylabel("abundance")
plt.plot(dat[:,0], dat[:,1], label="m=1 defensible")
plt.plot(dat[:,0], dat[:,2], label="m=1 efficient")
plt.plot(dat[:,0], dat[:,3], label="m=1 others")
plt.plot(dat[:,0], dat[:,4], label="CAPRI-3")
plt.plot(dat[:,0], dat[:,5], label="AON3")
plt.plot(dat[:,0], dat[:,6], label="m=3 FUSS")
# plt.legend()
plt.savefig(f"{data_dir}/ALL_S_abundance_30.pdf")

dat = np.loadtxt(f"{data_dir}/ALL_N_dep_beta_3.dat")

plt.ylim(-0.02, 1.02)
plt.xlabel("N")
plt.yticks(color="None")
#plt.ylabel("Abundance")
plt.plot(dat[:,0], dat[:,1], label="m=1 defensible")
plt.plot(dat[:,0], dat[:,2], label="m=1 efficient")
plt.plot(dat[:,0], dat[:,3], label="m=1 others")
plt.plot(dat[:,0], dat[:,4], label="CAPRI-3")
plt.plot(dat[:,0], dat[:,5], label="AON3")
plt.plot(dat[:,0], dat[:,6], label="m=3 FUSS")
plt.savefig(f"{data_dir}/ALL_N_dep_beta_3.pdf")

dat = np.loadtxt(f"{data_dir}/ALL_coop_level.dat")

plt.ylim(-0.02, 1.02)
plt.xlabel("benefit-to-cost ratio")
plt.ylabel("cooperation level")
plt.plot(dat[:,0], dat[:,28])
plt.savefig(f"{data_dir}/ALL_coop_level.pdf")


