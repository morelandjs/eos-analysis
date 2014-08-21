#!/usr/bin/env python
import numpy as np
import scipy 
import matplotlib.pyplot as plt
import prettyplotlib as ppl
from prettyplotlib import brewer2mpl

# plot WB EOS
pt, v2, v2err, v4, v4err = np.loadtxt("results/v2/WB-pions-20-30.dat",dtype='float',unpack=True)
#plt.errorbar(pt,v2,fmt="^-",yerr=v2err,color="blue",label="BW:")
plt.errorbar(pt,v4,fmt="^--",yerr=v4err,color="blue")

# plot HotQCD EOS
pt, v2, v2err, v4, v4err = np.loadtxt("results/v2/HotQCD-pions-20-30.dat",dtype='float',unpack=True)
#plt.errorbar(pt,v2,fmt="s-",yerr=v2err,color="red",label="HotQCD:")
plt.errorbar(pt,v4,fmt="s--",yerr=v4err,color="red")

# plot s95-PCE EOS
pt, v2, v2err, v4, v4err = np.loadtxt("results/v2/s95-PCE-pions-20-30.dat",dtype='float',unpack=True)
#plt.errorbar(pt,v2,fmt="o-",yerr=v2err,color="green",label="s95-PCE:")
plt.errorbar(pt,v4,fmt="o--",yerr=v4err,color="green")

# plot WB EOS
pt, v2, v2err, v4, v4err = np.loadtxt("results/v3/WB-pions-20-30.dat",dtype='float',unpack=True)
#plt.errorbar(pt,v2,fmt="^-",yerr=v2err,color="blue")
#plt.errorbar(pt,v4,fmt="^--",yerr=v4err,color="blue")

# plot HotQCD EOS
pt, v2, v2err, v4, v4err = np.loadtxt("results/v3/HotQCD-pions-20-30.dat",dtype='float',unpack=True)
#plt.errorbar(pt,v2,fmt="s-",yerr=v2err,color="red")
#plt.errorbar(pt,v4,fmt="s--",yerr=v4err,color="red")

# plot s95-PCE EOS
pt, v2, v2err, v4, v4err = np.loadtxt("results/v3/s95-PCE-pions-20-30.dat",dtype='float',unpack=True)
#plt.errorbar(pt,v2,fmt="o-",yerr=v2err,color="green")
#plt.errorbar(pt,v4,fmt="o--",yerr=v4err,color="green")

plt.xlabel("$p_T \,[GeV]$",fontsize=20)
plt.ylabel("$v_{2\{m\}}$",fontsize=25)
plt.annotate('$v\'_{2\{2\}}$', xy=(0.89, 0.87), xycoords='axes fraction',fontsize=20)
#plt.annotate('$v\'_{2\{4\}}$', xy=(0.9, 0.68), xycoords='axes fraction',fontsize=20)
plt.annotate('$v\'_{3\{2\}}$', xy=(0.89, 0.45), xycoords='axes fraction',fontsize=20)
#plt.annotate('$v\'_{3\{4\}}$', xy=(0.9, 0.31), xycoords='axes fraction',fontsize=20)
plt.annotate('pions: $20-30\%$', xy=(0.55, 0.05), xycoords='axes fraction',fontsize=25,color='gray')
plt.legend(loc='upper left', fontsize=17, frameon=False)
plt.xlim([0,2.75])
plt.ylim([0,0.23])
plt.show()
