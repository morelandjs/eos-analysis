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

fig = plt.figure(figsize=(columnwidth, columnwidth*aspect))
ax = fig.gca()

blue = plt.cm.Blues(0.6)
green = plt.cm.Greens(0.4)
red = plt.cm.Reds(0.5)

# hbarc
hbarc = 0.19733

# plot s95-PCE EOS
eEOS1, pEOS1, sEOS1, TEOS1 = np.loadtxt("s95-PCE-EOS.dat",dtype='float',unpack=True)

# plot WB EOS
eEOS2, pEOS2, sEOS2, TEOS2 = np.loadtxt("WB-EOS.dat",dtype='float',unpack=True)

# plot HotQCD EOS
eEOS3, pEOS3, sEOS3, TEOS3 = np.loadtxt("HotQCD-EOS.dat",dtype='float',unpack=True)

# plot HotQCD EOS from table (to check)
#T, I, p, e, s, Cv, cs = np.loadtxt("tables/EOS-table.dat",dtype='float',unpack=True, skiprows=1)

# compare energy density
plt.plot(TEOS1,eEOS1/TEOS1**4*hbarc**3, dashes = (5,2), linewidth=0.75, color=blue, label="s95-PCE")
plt.plot(TEOS2,eEOS2/TEOS2**4*hbarc**3, linewidth=0.75, color=green, label="BW")
plt.plot(TEOS3,eEOS3/TEOS3**4*hbarc**3, dashes = (20,2,5,2), linewidth=0.75, color=red,label="HotQCD")
#plt.plot(T*0.001,e,linewidth=1.0,color='purple',linestyle='--',label="HotQCD tbl")
#plt.plot(TEOS1,eEOS1,'--',linewidth=1.0,color='teal')
#plt.plot(TEOS2,eEOS2,linewidth=1.0,color='orange')

# compare entropy density
plt.plot(TEOS1,sEOS1/TEOS1**3*hbarc**3, dashes = (5,2), linewidth=0.75, color=blue)
plt.plot(TEOS2,sEOS2/TEOS2**3*hbarc**3, linewidth=0.75, color=green)
plt.plot(TEOS3,sEOS3/TEOS3**3*hbarc**3, dashes = (20,2,5,2), linewidth=0.75, color=red)
#plt.plot(T*0.001,s,linewidth=1.0,color='purple',linestyle='--',label="HotQCD tbl")
#plt.plot(TEOS1,sEOS1,'--',linewidth=1.0,color='teal')
#plt.plot(TEOS2,sEOS2,linewidth=1.0,color='orange')

# compare pressure
plt.plot(TEOS1,pEOS1/TEOS1**4*hbarc**3, dashes = (5,2), linewidth=0.75, color=blue)
plt.plot(TEOS2,pEOS2/TEOS2**4*hbarc**3, linewidth=0.75, color=green)
plt.plot(TEOS3,pEOS3/TEOS3**4*hbarc**3, dashes = (20,2,5,2), linewidth=0.75, color=red)
#plt.plot(T*0.001,p,linewidth=1.0,color='purple',linestyle='--',label="HotQCD tbl")
#plt.plot(TEOS1,pEOS1,'--',linewidth=1.0,color='teal')
#plt.plot(TEOS2,pEOS2,linewidth=1.0,color='orange')

#plt.plot(eEOS2,TEOS2,linewidth=1.0,color='orange')
#plt.plot(eEOS2,sEOS2,linewidth=1.0,color='orange')
#plt.plot(eEOS2,pEOS2,linewidth=1.0,color='orange')

#plt.legend(loc='upper right')
plt.xlim(0.0, 0.75)
plt.ylim(0.0, 24.0)
plt.xlabel('$T$ [GeV]', fontsize=texnormal)
plt.ylabel('')
#plt.title('Trace Anomaly')
#plt.savefig("trace.pdf")
#plt.yscale('log')
plt.legend(bbox_to_anchor=(1.02, 1.05), ncol=3, fontsize=texnormal)
plt.annotate('$s/T^3$', xy=(0.9, 0.775), xycoords='axes fraction',fontsize=texnormal)
plt.annotate('$e/T^4$', xy=(0.9, 0.575), xycoords='axes fraction',fontsize=texnormal)
plt.annotate('$p/T^4$', xy=(0.9, 0.175), xycoords='axes fraction',fontsize=texnormal)
#plt.annotate('dashed: VISHNU', xy=(0.45, 0.675), xycoords='axes fraction',fontsize=15,color='gray')
#plt.annotate('solid orange: WB', xy=(0.45, 0.575), xycoords='axes fraction',fontsize=15,color='gray')
#plt.annotate('solid red: WB', xy=(0.45, 0.475), xycoords='axes fraction',fontsize=15,color='gray')
despine()
plt.tight_layout(pad=0, w_pad=.1)
plt.savefig("eos_comparison.pdf")
