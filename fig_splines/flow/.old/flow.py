#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import glob, os, sys

# loop over equations of state
for eos, path in zip([os.path.basename(path) for path in sys.argv[1:]], sys.argv[1:]):
       
    print("parsing", eos + '/20to30')
    filenames = glob.glob(path + '/20to30/*pi.dflow.gz')
    pt, v2_, v3_ = np.array([np.loadtxt(f, usecols=(0,8,11)) for f in filenames]).transpose(2,0,1)
    pt = pt[0]
    v2 = np.mean(v2_, axis=0)
    dv2 = np.std(v2_, axis=0)/np.sqrt(v2_.size)
    v3 = np.mean(v3_, axis=0)
    dv3 = np.std(v3_, axis=0)/np.sqrt(v3_.size)

    directory = "results/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    np.savez_compressed(directory + eos + ".npz", pt=pt, v2=v2, dv2=dv2, v3=v3, dv3=dv3)
