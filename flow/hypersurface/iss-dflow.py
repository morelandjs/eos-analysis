#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import glob, os, sys

# specify EOS folder, e.g. path/HotQCD-errors/, path/HotQCD, path/s95-PCE, etc
for centrality in "/0to10/", "/20to30/", "/40to50/":
        
    eos = os.path.basename(sys.argv[1])
    particle_types = "/*.pi.dflow.gz", "/*.k.dflow.gz", "/*.p.dflow.gz"
    names = "/pion", "/kaon", "proton"
    print("parsing", eos + centrality)
    
    for particle_type, name in zip(particle_types, names):
        filenames = glob.glob(sys.argv[1] + centrality + particle_type)
        pt, v2_, v3_ = np.array([np.loadtxt(f, usecols=(0,8,11)) for f in filenames]).transpose(2,0,1)
        pt = pt[0]
        v2 = np.mean(v2_, axis=0)
        dv2 = np.std(v2_, axis=0)/np.sqrt(v2_.size)
        v3 = np.mean(v3_, axis=0)
        dv3 = np.std(v3_, axis=0)/np.sqrt(v3_.size)

        directory = "../results/" + eos + centrality
        if not os.path.exists(directory):
            os.makedirs(directory)
        np.savez_compressed(directory + name, pt=pt, v2=v2, dv2=dv2, v3=v3, dv3=dv3)
