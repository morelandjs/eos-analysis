#!/usr/bin/env python

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy 

sns.set(font='CMU Serif', font_scale=1.5, rc={"lines.linewidth": 2.5})

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
plt.plot(TEOS1,eEOS1/TEOS1**4*hbarc**3,'--',linewidth=1.0,color='teal',label="s95-PCE")
plt.plot(TEOS2,eEOS2/TEOS2**4*hbarc**3,linewidth=1.0,color='orange',label="BW")
plt.plot(TEOS3,eEOS3/TEOS3**4*hbarc**3,linewidth=1.0,color='red',label="HotQCD")
#plt.plot(T*0.001,e,linewidth=1.0,color='purple',linestyle='--',label="HotQCD tbl")
#plt.plot(TEOS1,eEOS1,'--',linewidth=1.0,color='teal')
#plt.plot(TEOS2,eEOS2,linewidth=1.0,color='orange')

# compare entropy density
plt.plot(TEOS1,sEOS1/TEOS1**3*hbarc**3,'--',linewidth=1.0,color='teal')
plt.plot(TEOS2,sEOS2/TEOS2**3*hbarc**3,linewidth=1.0,color='orange')
plt.plot(TEOS3,sEOS3/TEOS3**3*hbarc**3,linewidth=1.0,color='red')
#plt.plot(T*0.001,s,linewidth=1.0,color='purple',linestyle='--',label="HotQCD tbl")
#plt.plot(TEOS1,sEOS1,'--',linewidth=1.0,color='teal')
#plt.plot(TEOS2,sEOS2,linewidth=1.0,color='orange')

# compare pressure
plt.plot(TEOS1,pEOS1/TEOS1**4*hbarc**3,'--',linewidth=1.0,color='teal')
plt.plot(TEOS2,pEOS2/TEOS2**4*hbarc**3,linewidth=1.0,color='orange')
plt.plot(TEOS3,pEOS3/TEOS3**4*hbarc**3,linewidth=1.0,color='red')
#plt.plot(T*0.001,p,linewidth=1.0,color='purple',linestyle='--',label="HotQCD tbl")
#plt.plot(TEOS1,pEOS1,'--',linewidth=1.0,color='teal')
#plt.plot(TEOS2,pEOS2,linewidth=1.0,color='orange')

#plt.plot(eEOS2,TEOS2,linewidth=1.0,color='orange')
#plt.plot(eEOS2,sEOS2,linewidth=1.0,color='orange')
#plt.plot(eEOS2,pEOS2,linewidth=1.0,color='orange')

#plt.legend(loc='upper right')
plt.xlim(0.0, 0.75)
plt.ylim(0.0, 24.0)
plt.xlabel('$T$ [MeV]')
plt.ylabel('')
#plt.title('Trace Anomaly')
#plt.savefig("trace.pdf")
#plt.yscale('log')
plt.legend(loc=9,ncol=3,fontsize=15)
plt.annotate('$s/T^3$', xy=(0.9, 0.775), xycoords='axes fraction',fontsize=15)
plt.annotate('$e/T^4$', xy=(0.9, 0.575), xycoords='axes fraction',fontsize=15)
plt.annotate('$p/T^4$', xy=(0.9, 0.175), xycoords='axes fraction',fontsize=15)

#plt.annotate('dashed: VISHNU', xy=(0.45, 0.675), xycoords='axes fraction',fontsize=15,color='gray')
#plt.annotate('solid orange: WB', xy=(0.45, 0.575), xycoords='axes fraction',fontsize=15,color='gray')
#plt.annotate('solid red: WB', xy=(0.45, 0.475), xycoords='axes fraction',fontsize=15,color='gray')
plt.show()
