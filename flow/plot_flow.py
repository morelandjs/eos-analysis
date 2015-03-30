#!/usr/bin/env python
import numpy as np
import scipy 
import matplotlib.pyplot as plt
import sys

for f in sys.argv[1:]:
    pt = np.arange(0.2,2.8,0.2) + 0.1
    vn2, vn2_err, vn4, vn4_err = np.loadtxt(f,skiprows=1,unpack=True)
    plt.errorbar(pt,vn2,yerr=vn2_err)
    plt.errorbar(pt,vn4,yerr=vn4_err)
plt.show()
