#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import math
import sys

# plot s95-PCE spectra
###################################
pt = np.loadtxt("spectra-files/00-05/s95-kaon-pt.dat",unpack=True)
y,binEdges=np.histogram(pt,100,weights=1/(4.*np.pi*pt),normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.plot(bincenters,y,'-', label="s95-PCE")
#plt.hist(pt,100, weights=1/(4.*np.pi*pt), normed=True, histtype='step', label="s95-PCE")

# plot WB spectra
###################################
pt = np.loadtxt("spectra-files/00-05/WB-kaon-pt.dat",unpack=True)
y,binEdges=np.histogram(pt,100,weights=1/(4.*np.pi*pt),normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.plot(bincenters,y,'-', label="WB")
#plt.hist(pt,100, weights=1/(4.*np.pi*pt), normed=True, histtype='step', label="WB")

# plot Hot-QCD spectra
###################################
pt = np.loadtxt("spectra-files/00-05/HotQCD-kaon-pt.dat",unpack=True)
y,binEdges=np.histogram(pt,100,weights=1/(4.*np.pi*pt),normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.plot(bincenters,y,'-', label="HotQCD")
#plt.hist(pt,100, weights=1/(4.*np.pi*pt), normed=True, histtype='step', label="HotQCD")

# plot properties
###################################################################
plt.yscale('log')
plt.xlim([0,5.0])
plt.ylim([10e-5,1])
plt.xlabel('$p_T$ [GeV]',fontsize=15)
plt.ylabel(r'$\frac{1}{2 \pi p_T} \, \frac{dN}{dy dp_T}$',fontsize=20)
plt.legend(loc='lower left',fontsize=15,frameon=False)
plt.annotate('$\pi$', xy=(0.9, 0.9), xycoords='axes fraction',fontsize=40,color='black')
plt.show()
