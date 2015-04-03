#!/usr/bin/env python
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy
import sys

sns.set(font='CMU Serif', font_scale=1.5, rc={"lines.linewidth": 2.5})

for f in sys.argv[1:]:
    # plot UrQMD data
    #pt = np.arange(0.2,2.8,0.2) + 0.1
    #vn2, vn2_err, vn4, vn4_err = np.loadtxt(f,skiprows=1,unpack=True)
    #plt.plot(pt,vn2,linewidth=1)
    #plt.errorbar(pt,vn2,yerr=vn2_err,linewidth=1)
    #plt.errorbar(pt,vn4,yerr=vn4_err)

    # plot iSS data
    pt, v2, dv2, v3, dv3 = np.loadtxt(f,unpack=True)
    plt.plot(pt, v2, linewidth=1)
    plt.plot(pt, v3, linewidth=1)
    #plt.fill_between(pt, v2 - dv2, v2 + dv2, alpha=0.25)
    #plt.fill_between(pt, v3 - dv3, v3 + dv3, alpha=0.25)

#plt.annotate('$v_2\{2\}$', xy=(0.8,0.9), xycoords='axes fraction', fontsize=16, ha='center')
#plt.annotate('$v_3\{2\}$', xy=(0.8,0.55), xycoords='axes fraction', fontsize=16, ha='center')
plt.annotate('$v_2$', xy=(0.8,0.73), xycoords='axes fraction', fontsize=20, ha='center')
plt.annotate('$v_3$', xy=(0.8,0.35), xycoords='axes fraction', fontsize=20, ha='center')
plt.annotate('$\pi^+$', xy=(0.2,0.88), xycoords='axes fraction', fontsize=40, ha='center')
#plt.annotate('$0.2 < p_T < 3.0$, $|\eta|<2.5$', xy=(0.2,0.82), xycoords='axes fraction', fontsize=18, ha='center')
plt.annotate('$20-30\%$', xy=(0.8,0.05), xycoords='axes fraction', fontsize=30, ha='center')
plt.xlabel(r"$p_T$ [GeV]")
#plt.ylabel(r"$v_n\{2\}$")
plt.ylabel(r"$v_n$")
plt.xlim([0,3.0])
plt.ylim([0,0.25])
plt.show()
