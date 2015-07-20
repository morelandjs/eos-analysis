#!/usr/bin/env python3

from __future__ import division, print_function, unicode_literals

import functools

import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib import cm
import numpy as np
import scipy
import seaborn as sns
import os, sys

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

fig = plt.figure(figsize=(1.3*columnwidth, columnwidth*aspect))
ax = fig.gca()

blue = '#3399FF'
magenta = '#B8008A'
orange = '#FF9900'

# hbarc
hbarc = 0.19733

# loa s95, HotQCD and WB equations of state
eEOS1, pEOS1, sEOS1, TEOS1 = np.loadtxt("s95-PCE-EOS.dat",dtype='float',unpack=True)
eEOS2, pEOS2, sEOS2, TEOS2 = np.loadtxt("HotQCD-EOS.dat",dtype='float',unpack=True)
eEOS3, pEOS3, sEOS3, TEOS3 = np.loadtxt("WB-EOS.dat",dtype='float',unpack=True)

# compare energy density
plt.plot(TEOS1, eEOS1/TEOS1**4*hbarc**3, linewidth=0.75, color=blue, label="S95", zorder=3)
plt.plot(TEOS2[0:155500], eEOS2[0:155500]/TEOS2[0:155500]**4*hbarc**3, linewidth=0.75, color=magenta, label="HotQCD", zorder=2)
plt.plot(TEOS3[0:155500], eEOS3[0:155500]/TEOS3[0:155500]**4*hbarc**3, linewidth=0.75, color=orange,label="WB", zorder=1)

# compare entropy density
plt.plot(TEOS1,sEOS1/TEOS1**3*hbarc**3, linewidth=0.75, color=blue, zorder=3)
plt.plot(TEOS2[0:155500], sEOS2[0:155500]/TEOS2[0:155500]**3*hbarc**3, linewidth=0.75, color=magenta, zorder=2)
plt.plot(TEOS3[0:155500], sEOS3[0:155500]/TEOS3[0:155500]**3*hbarc**3, linewidth=0.75, color=orange, zorder=1)

# compare pressure
plt.plot(TEOS1,pEOS1/TEOS1**4*hbarc**3, linewidth=0.75, color=blue, zorder=3)
plt.plot(TEOS2[0:155500], pEOS2[0:155500]/TEOS2[0:155500]**4*hbarc**3, linewidth=0.75, color=magenta, zorder=2)
plt.plot(TEOS3[0:155500], pEOS3[0:155500]/TEOS3[0:155500]**4*hbarc**3, linewidth=0.75, color=orange, zorder=1)


plt.xlim(0.0, 0.75)
plt.ylim(0.0, 24.0)
plt.xlabel('$T$ [GeV]', fontsize=texnormal)
plt.ylabel('')
plt.legend(bbox_to_anchor=(0.91, 1.05), ncol=3, handlelength=2.9, fontsize=texnormal)
plt.annotate('$s/T^3$', xy=(0.9, 0.775), xycoords='axes fraction',fontsize=texnormal)
plt.annotate('$e/T^4$', xy=(0.9, 0.575), xycoords='axes fraction',fontsize=texnormal)
plt.annotate('$p/T^4$', xy=(0.9, 0.175), xycoords='axes fraction',fontsize=texnormal)
despine()
plt.tight_layout(pad=0, w_pad=.1)
plt.savefig("eos_comparison.pdf")
