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
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['axes.labelsize'] = 30
plt.rcParams['font.size'] = 22
plt.rcParams['axes.linewidth'] = 1.2

# +
# take histogram
dat = np.loadtxt("n4_samples.dat")

benefit = 3.0
p_range = [benefit/4-1.0, 3*benefit/4]

Hb, xedges, yedges = np.histogram2d(dat[:,0], dat[:,1], bins=100, range=[p_range, p_range])
Hc, xedges, yedges = np.histogram2d(dat[:,0], dat[:,2], bins=100, range=[p_range, p_range])
Hd, xedges, yedges = np.histogram2d(dat[:,0], dat[:,3], bins=100, range=[p_range, p_range])
H4 = (Hb + Hc + Hd).T # Let each row list bins with common y range.


# +
# feasible set
# dddd : 0, 0
# dddc : b/4, [b/4-1,b/4]
# ddcc : 2b/4, [2b/4-1, 2b/4]
# dccc : 3b/4, 3b/4-1
# cccc : b-1, b-1
# cccd : 3b/4-1, [3b/4-1, 3b/4]
# ccdd : 2b/4-1, [2b/4-1, 2b/4]
# cddd :  b/4-1, b/4

def make_fig(histo, b):
    c = 'gray'
    lw = 0.5
    p0 = (0, 0)
    p1 = (b/4, b/4-1)
    p2 = (2*b/4, 2*b/4-1)
    p3 = (3*b/4, 3*b/4-1)
    p2x= (2*b/4, 2*b/4)
    p4 = (b-1, b-1)
    p5 = (3*b/4-1, 3*b/4)
    p6 = (2*b/4-1, 2*b/4)
    p7 = (b/4-1, b/4)
    p_range = (b/4-1, 3*b/4)
    def draw_line(b, e, ls='solid'):
        x = np.linspace(b[0], e[0], 100)
        y = np.linspace(b[1], e[1], 100)
        plt.plot(x,y,color=c, linewidth=0.5, linestyle=ls)
    draw_line(p0, p1)
    draw_line(p1, p2)
    draw_line(p2, p3)
    draw_line(p3, p4)
    draw_line(p4, p5)
    draw_line(p5, p6)
    draw_line(p6, p7)
    draw_line(p7, p0)
    draw_line(p3, p2x, ls='dotted')
    draw_line(p0, p4, ls='dotted')
    plt.xlim(p_range[0], p_range[1])
    plt.ylim(p_range[0], p_range[1])
    plt.plot(2, 2, 'o', color='magenta')
    plt.plot(0, 0, 'x', color='black')
    plt.imshow(histo, cmap='Blues', origin='lower', interpolation='antialiased', extent=(p_range[0],p_range[1],p_range[0],p_range[1]), norm=colors.LogNorm())

make_fig(H4,benefit)
plt.savefig("n4_dist.pdf")

# +
# when n=3, benefit = 1.5
# payoff range : [-0.5, 1.0]
dat = np.loadtxt("n3_samples.dat")
benefit = 2.0
p_range = [benefit * 1/3 - 1, benefit * 2/3]
Hb, xedges, yedges = np.histogram2d(dat[:,0], dat[:,1], bins=100, range=[p_range,p_range])
Hc, xedges, yedges = np.histogram2d(dat[:,0], dat[:,2], bins=100, range=[p_range,p_range])
H3 = (Hb + Hc).T # Let each row list bins with common y range.


# -

def draw_theoretical_n3(b):
    c = 'gray'
    p0 = (0.0, 0.0)
    p1 = (b * 1/3, b * 1/3 - 1)
    p1x= (b * 1/3, b * 1/3)
    p2 = (b * 2/3, b * 2/3 -1)
    p3 = (b - 1, b - 1)
    p4 = (b * 2/3 - 1, b * 2/3)
    p5 = (b * 1/3 - 1, b * 1/3)
    def draw_line(s, e, ls='solid'):
        x = np.linspace(s[0], e[0], 100)
        y = np.linspace(s[1], e[1], 100)
        plt.plot(x,y,color=c, linewidth=0.5, linestyle=ls)
    draw_line(p0, p1)
    draw_line(p1, p2)
    draw_line(p2, p3)
    draw_line(p3, p4)
    draw_line(p4, p5)
    draw_line(p5, p0)
    draw_line(p0, p3, ls='dotted')
    draw_line(p2, p1x, ls='dotted')
    plt.xlim(p5[0], p2[0])
    plt.ylim(p1[1], p4[1])
    plt.xticks([0,1])
    plt.yticks([0,1])
    plt.plot(1, 1, 'o', color='magenta')
    plt.plot(0, 0, 'x', color='black')
draw_theoretical_n3(benefit)
plt.imshow(H3, cmap='Blues', origin='lower', interpolation='antialiased', extent=(p_range[0],p_range[1],p_range[0],p_range[1]), norm=colors.LogNorm())
plt.savefig("n3_dist.pdf")

# when n=2
# payoff range : [-0.5, 1.0]
dat = np.loadtxt("n2_capri2_samples.dat")
benefit = 1.5
p_range = [benefit/2 - 1, benefit/2]
Hb, xedges, yedges = np.histogram2d(dat[:,0], dat[:,1], bins=100, range=[p_range,p_range])
H2 = Hb.T # Let each row list bins with common y range.


def draw_theoretical_n2(b):
    c = 'gray'
    p0 = (0.0, 0.0)
    p1 = (b * 1/2, b * 1/2 - 1)
    p2 = (b - 1, b - 1)
    p3 = (b * 1/2 - 1, b * 1/2)
    def draw_line(s, e, ls='solid'):
        x = np.linspace(s[0], e[0], 100)
        y = np.linspace(s[1], e[1], 100)
        plt.plot(x,y,color=c, linewidth=0.5, linestyle=ls)
    draw_line(p0, p1)
    draw_line(p1, p2)
    draw_line(p2, p3)
    draw_line(p3, p0)
    draw_line(p0, p2, ls='dotted')
    plt.xlim(p3[0], p1[0])
    plt.ylim(p1[1], p3[1])
    plt.xticks([0,0.5])
    plt.yticks([0,0.5])
    plt.plot(p2[0], p2[1], 'o', color='magenta')
    plt.plot((p0[0]+p2[0])/2, (p0[1]+p2[1])/2, 'x', color='black')
draw_theoretical_n2(benefit)
plt.imshow(H2, cmap='Blues', origin='lower', interpolation='antialiased', extent=(p_range[0],p_range[1],p_range[0],p_range[1]), norm=colors.LogNorm())
plt.savefig("n2_dist.pdf")


