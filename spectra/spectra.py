#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import math
import sys

# scientific dictionary
scientific = {i: lambda x: float(x.replace('D', 'E')) for i in range(0,23)}

species_min = 101 # pion
species_max = 101 # pion
#species_min = 1 # proton 
#species_max = 1 # proton 
#species_min = 27 # lambda
#species_max = 39 # lambda
#species_min = 106 # kaon
#species_max = 106 # kaon

# open output file
with open('HotQCD-pion-pt.dat','w') as wf:

    # event number
    ievent = 0

    # parse spectra (avg over events)
    for f in sys.argv[1:]:
        E, px, py, pz, m, ityp = np.loadtxt(f,usecols=(4,5,6,7,8,9),skiprows=16,unpack=True,converters=scientific)
        E = np.asarray(E)
        px = np.asarray(px)
        py = np.asarray(py)
        pz = np.asarray(py)
        y = np.log((E+pz)/(E-pz))/2.
        ityp = np.asarray(ityp)
        N = ityp.size
        px = np.asarray([px[i] for i in range(N) if ityp[i] >= species_min and ityp[i] <= species_max and np.abs(y[i]) < 1.0])
        py = np.asarray([py[i] for i in range(N) if ityp[i] >= species_min and ityp[i] <= species_max and np.abs(y[i]) < 1.0])
        pz = np.asarray([pz[i] for i in range(N) if ityp[i] >= species_min and ityp[i] <= species_max and np.abs(y[i]) < 1.0])
        ityp = np.asarray([ityp[i] for i in range(N) if ityp[i] >= species_min and ityp[i] <= species_max and np.abs(y[i]) < 1.0])
        pt = np.sqrt(px**2.+py**2.) 
        for item in np.nditer(pt):
            wf.write("%s\n" % item)
        print ievent
        ievent += 1
