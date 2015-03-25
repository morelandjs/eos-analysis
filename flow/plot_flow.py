#!/usr/bin/env python
import numpy as np
import scipy 
import matplotlib.pyplot as plt
import sys

for f in sys.argv[1:]:
    print(f)
    pt = np.arange(0.0,2.0,0.25) + 0.125
    vn2, vn4 = np.loadtxt(f,skiprows=1,unpack=True)
    plt.plot(pt,vn2)
    plt.plot(pt,vn4)
    #plt.errorbar(pt, vn2, yerr=vn2err)
    #plt.errorbar(pt, vn4, yerr=vn4err)
plt.show()
