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

data_dir = "evo_fixation_probs"

dat1 = np.loadtxt(f"{data_dir}/N10_b2_n2_random.dat")
dat2 = np.loadtxt(f"{data_dir}/N10_b2_n2_capri2.dat")

plt.yscale('log')
plt.ylim(1.0e-1, 1.0e3)
plt.xlim([0,0.3])
plt.xlabel(r"$\rho$")
plt.ylabel(r'$P(\rho)$')
plt.plot(dat1[:,0]+0.001, dat1[:,1]*1000, label="Random")
plt.plot(dat2[:,0]+0.001, dat2[:,1]*1000, label="CAPRI-2")
plt.legend()
plt.plot([0.1,0.1],[1.0e-1,1.0e3],'--',color='gray')
plt.savefig("rho_distribution_n2.pdf")

dat3 = np.loadtxt(f"{data_dir}/N10_b2_n3_random.dat")
dat4 = np.loadtxt(f"{data_dir}/N10_b2_n3_capri3.dat")

plt.yscale('log')
plt.ylim(1.0e-1, 1.0e3)
plt.yticks(color="None")
plt.xlim([0,0.3])
plt.xlabel(r"$\rho$")
plt.plot(dat3[:,0]+0.001, dat3[:,1]*1000, label="Random")
plt.plot(dat4[:,0]+0.001, dat4[:,1]*1000, label="CAPRI-3", linewidth=2.0)
plt.legend()
plt.plot([0.1,0.1],[1.0e-1,1.0e3],'--',color='gray')
plt.savefig("rho_distribution_n3.pdf")


