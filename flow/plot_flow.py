#!/usr/bin/env python
import numpy as np
import scipy 
import matplotlib.pyplot as plt
import sys

for f in sys.argv[1:]:
    print(f)
    pt = np.arange(0.125,2.625,0.25)
    print(pt.size)
    vn2, vn2err, vn4, vn4err = np.loadtxt(f,skiprows=3,unpack=True)
    plt.errorbar(pt, vn2, yerr=vn2err)
    #plt.errorbar(pt, vn4, yerr=vn4err)
plt.show()
