#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
import glob

# figure properties
aspect = 1/1.618
resolution = 72.27
columnwidth = 246/resolution
textwidth = 510/resolution
textiny, texsmall, texnormal = 8., 9.25, 10.
offblack = '#262626'

plt.rcdefaults()
plt.rcParams.update({
        'font.family': ['CMU Serif'],
        'font.size': texsmall,
        'legend.fontsize': 8.5,
        'axes.labelsize': texsmall,
        'axes.titlesize': texsmall,
        'xtick.labelsize': textiny,
        'ytick.labelsize': textiny,
        'lines.linewidth': .9,
        'patch.linewidth': .9,
        'axes.linewidth': .5,
        'text.color': offblack,
        'axes.edgecolor': offblack,
        'axes.labelcolor': offblack,
        'xtick.color': offblack,
        'ytick.color': offblack,
        'xtick.major.size': 2,
        'ytick.major.size': 2,
        'xtick.major.pad': 1.8,
        'ytick.major.pad': 1.8,
        'legend.numpoints': 1,
        'legend.scatterpoints': 1,
        'legend.frameon': False,
        'pdf.fonttype': 42
})

def despine(ax=None, remove_ticks=False):
    if ax is None:
        ax = plt.gca()
    for spine in 'top', 'right':
        ax.spines[spine].set_visible(False)
    if remove_ticks:
        for ax_name in 'xaxis', 'yaxis':
            getattr(ax, ax_name).set_ticks_position('none')
    else:
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

blue = '#3399FF'
magenta = '#B8008A'
orange = '#FF9900'
green = '#4DB84D'
gray = '#504A4B'
hbarc = 0.19733

fig, axes = plt.subplots(2, figsize=(columnwidth, 1.6*columnwidth*aspect), sharex=True)

# HotQCD equation of state
e, p, s, T = np.loadtxt('hotqcd.dat', unpack=True)
hqcd = Spline(e, p)

# Wuppertal-Budapest equation of state
e, p, s, T = np.loadtxt('wb.dat', unpack=True)
wb = Spline(e, p)

# s95p-v1 equation of state
e, p, s, T = np.loadtxt('s95.dat', unpack=True)
s95 = Spline(e, p)

for ax, length, path, label in zip(axes, [0.8, 1], ['bayesian/prior/', 'bayesian/posterior/'], ['Prior', 'Posterior']):
    # HotQCD equation of state
    ax.plot(T, hqcd(e, nu=1), linewidth=1, color=plt.cm.Blues(0.6), label='HQCD matched', zorder=3)

    # Wuppertal-Budapest equation of state
    ax.plot(T, wb(e, nu=1), linewidth=1, color=plt.cm.Reds(0.6), label='WB matched', zorder=4)

    # s95p-v1 equation of state
    ax.plot(T, s95(e, nu=1), linewidth=1, color=plt.cm.Greens(0.6), label=r'S95 matched', zorder=2)

    # scott pratt prior/posterior
    ax.plot([0], np.zeros(1), linewidth=0.5, color=gray, zorder=1, label='Bayesian')
    for index, f in enumerate(glob.glob(path + '*.dat')):
        T_, cs2 = np.loadtxt(f, usecols=(0,3)).T
        ax.plot(T_/1000, cs2, linewidth=0.2, color='0.4', alpha=0.4, zorder=1)

    ax.annotate(label, xy=(0.25, 0.33), ha='center', va='bottom')
    
    # axis properties
    ax.set_xlim(0.1,0.4)
    ax.set_xticks([0.1, 0.2, 0.3, 0.4])
    ax.set_ylim(0, 0.38)
    ax.set_yticks([0, 0.1, 0.2, 0.3]) if ax == axes[1] else ax.set_yticks([0.1, 0.2, 0.3])
    ax.set_xlabel('$T$ [GeV]') if ax == axes[1] else ax.set_xlabel('')
    ax.set_ylabel('$c_s^2$')
    despine(ax)

# figure properties
#################################################################
plt.tight_layout(pad=0)
handles, labels = axes[0].get_legend_handles_labels()
axes[1].legend(bbox_to_anchor=(1.02,0.6))
plt.savefig('cs.pdf')
#plt.show()
