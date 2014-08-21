#!/usr/bin/env python

import matplotlib.pyplot as plt
import itertools
from itertools import izip
import numpy as np
import math
import sys
import csv

# scientific dictionary
scientific = {i: lambda x: float(x.replace('D', 'E')) for i in range(0,23)}

# event counter
nev = 0

# open output file
with open('tmp/HotQCD-dflow-20-30.dat','wb') as wf:
    
    # parse spectra (avg over events)
    for f in sys.argv[1:]:

        # load an event and store (ityp,phi,pt,y) into an array
        E, px, py, pz, ityp = np.loadtxt(f,usecols=(4,5,6,7,9),skiprows=16,unpack=True,converters=scientific)
        ityp = np.asarray(ityp)
        phi = np.arctan2(np.asarray(py),np.asarray(px))
        pt = np.sqrt(np.asarray(px)**2.+np.asarray(py)**2.)
        y = np.log((np.asarray(E)+np.asarray(pz))/(np.asarray(E)-np.asarray(pz)))/2.
        
        # cut the data on rapidity 
        data = np.vstack((ityp,phi,pt,y)).T 
        data = data[(data[:,3] > -1.0) & (data[:,3] <  1.0)] 

        # print output to file
        for row in data:
            wf.write("{:20}{:20}{:20}\n".format(*row))

        # skip line   
        print "event: ",nev
        wf.write("\n")
        nev += 1

