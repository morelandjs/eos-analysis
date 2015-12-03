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
gray = '#C9C9C9'

# UrQMD equation of state
#################################################################

T, I, e, p, Cs = np.loadtxt("tables/urqmd/hrg-urqmd-eos.dat",dtype='float',unpack=True,skiprows=1)
plt.plot(T/1000, I, dashes=(4,1), color=offblack,linewidth=1, zorder = 301, label="HRG")

# HotQCD equation of state
#################################################################

# temperature limits [GeV]
Tmin = 0.130
Tmax = 0.800
dT = 0.001
nT = int((Tmax-Tmin)/dT)
T = np.linspace(Tmin,Tmax,nT)

# eos parameters
hbarc = 0.19733
Tc = 0.154
ct = 3.8706
an = -8.7704
bn = 3.9200
cn = 0.0
dn = 0.3419
t0 = 0.9761
ad = -1.2600
bd = 0.8425
cd = 0.0
dd = -0.0475
tb = T/Tc
pid = 47.5*np.pi**2./90.

# pressure over T^4
p = 0.5*(1.+np.tanh(ct*(tb-t0)))*(pid+an/tb+bn/tb**2.+cn/tb**3.+dn/tb**4.)/(1.+ad/tb+bd/tb**2.+cd/tb**3.+dd/tb**4.)
dfp = interp1d(T, interpolate.splev(T,interpolate.splrep(T,p,s=0),der=1), kind='cubic')
fp = interp1d(T, p, kind='cubic')

# trace anomaly over T^4
I = []
for iT in T:
    I.append(dfp(iT)*iT)
fI_hi = interp1d(T*1000, I, kind='cubic')

# plot HotQCD trace anomaly))
plt.plot(T, I, color=plt.cm.Blues(0.6),linewidth=1, zorder = 99, label="HQCD")

# Wuppertal-Budapest equation of state
#################################################################

# temperature limits in GeV
T0 = 0.2
Tmin = 0.100
Tmax = 0.800
dT = 0.005
nT = int((Tmax-Tmin)/dT)
Tvec = np.linspace(Tmin,Tmax,nT)
t=Tvec/T0

# eos parameters
hbarc=0.19733
h0 = 0.1396
h1 = -0.1800
h2 = 0.0350 
f0 = 1.05 
f1 = 6.39 
f2 = -4.72 
g1 = -0.92 
g2 = 0.57

# trace anomaly: I/T^4
I = np.exp(-h1/t-h2/t**2.)*(h0 + f0*(np.tanh(f1*t+f2)+1.)/(1.+g1*t+g2*t**2.)) 
I = np.asarray(I)
plt.plot(Tvec, I, color=plt.cm.Reds(0.6), label='WB')

# s95p-v1 equation of state
#################################################################

T = np.loadtxt("tables/s95p-v1/s95p-v1_Teospar.dat", usecols=(0,), skiprows=1, unpack=True)/1000.
e, p = np.loadtxt("tables/s95p-v1/s95p-v1_Teosdens.dat", usecols=(0,1,), skiprows=1, unpack=True)*hbarc**3
plt.plot(T, (e-3.*p)/T**4, linewidth=1, color=plt.cm.Greens(0.6), label='s95p-v1')

# figure properties
#################################################################
despine()
plt.xlim(0,0.63)
plt.ylim(0,6)
plt.yticks([1,2,3,4,5,6])
plt.xlabel('$T$ [GeV]')
plt.ylabel('$(e-3p)/T^4$')
plt.tight_layout(pad=0.2)
plt.legend(bbox_to_anchor=(1,1.05))
plt.savefig('trace.pdf')
plt.show()
