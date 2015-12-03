#!/usr/bin/env python3

from __future__ import division, print_function, unicode_literals

import functools

import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib import cm
import numpy as np
import scipy
from scipy import interpolate
from scipy.interpolate import interp1d
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

blue = '#3399FF'
magenta = '#B8008A'
orange = '#FF9900'
green = '#4DB84D'
gray = '#5C5858'

hbarc = 0.19733

####################################################################################

# UrQMD interaction measure 

T, I, e, p, Cs = np.loadtxt("hrg-eos/hrg-urqmd-eos.dat",dtype='float',unpack=True,skiprows=1)
m = len(T)
UrQMD = interp1d(T/1000., I, kind='cubic')
plt.plot(T/1000., I, linewidth=1, color=offblack, dashes=(4,1), zorder=3, label='HRG')

# HotQCD equation of state
#################################################################

e, p, s, T = np.loadtxt('eos-tables/hotqcd.dat', unpack=True)
e *= hbarc**3
p *= hbarc**3
plt.plot(T, (e-3.*p)/T**4, linewidth=1, color=plt.cm.Blues(0.8), dashes=(8,2), zorder=3, label=r'HQCD matched')

# HotQCD splines
#################################################################
plt.plot(np.ones(10), -np.ones(10), linewidth=0.5, color=plt.cm.bone(0.5), alpha=1, label=r'HQCD samples', zorder=1)
for index, f in enumerate('eos-tables/hotqcd_{}.dat'.format(i) for i in range(100)):
    e, p, s, T = np.loadtxt(f, unpack=True)
    e *= hbarc**3
    p *= hbarc**3
    plt.plot(T, (e-3.*p)/T**4, linewidth=0.2, color=plt.cm.bone(index/100), alpha=0.3, zorder=1)

# hydro-to-urqmd transition line
#################################################################
plt.axvline(x=0.155, ymin=0., ymax = 0.9, linewidth=0.2, color=offblack)
plt.annotate(r'$0.155$ GeV', xy=(0.155, 5.6), va='bottom', ha='center', fontsize=textiny)

# figure properties
#################################################################
despine()
plt.xlim(0,0.45)
plt.xticks([0, 0.1, 0.2, 0.3, 0.4])
plt.ylim(0,6)
plt.yticks([1,2,3,4,5,6])
plt.xlabel('$T$ [GeV]')
plt.ylabel('$(e-3p)/T^4$')
plt.tight_layout(pad=0.2)
plt.legend(bbox_to_anchor=(1,1.05))
plt.savefig('splines.pdf')
#plt.show()
