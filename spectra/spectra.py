#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import itertools
import subprocess
import sys
import glob, os

# instructions run on list of folder directories, e.g. ./spectra.py HotQCD-EOS1 HotQCD-EOS2 ... HotQCD-EOS99

# event generator
def events(filename):
    event = []
    with subprocess.Popen(('zcat', filename), stdout=subprocess.PIPE) as proc:
        for l in proc.stdout:
            if len(l) == 435:
                p = l[97:192].replace(b'D', b'E')
                E = p[0:23]
                px = p[24:47]
                py = p[48:71]
                pz = p[72:95]
                ityp = l[218:221]
                chg = l[225:227]
                event.append((E, px, py, pz, ityp, chg))
            else:
                if event:
                    yield np.array(event, dtype=float)
                    event.clear()

        if event:
            yield np.array(event, dtype=float)


# specify EOS folder, e.g. path/HotQCD-errors/, path/HotQCD, path/s95-PCE, etc
for centrality in "/0to10/", "/20to30/", "/40to50/":
    
    eos = os.path.basename(sys.argv[1])
    filenames = glob.glob(sys.argv[1] + centrality + "/*.urqmd.gz")
    print(centrality)

    pion = []
    kaon = []
    proton = []
    nev = 0

    # loop over events in each folder
    for event in itertools.chain.from_iterable(events(f) for f in filenames[:]):

        # unpack particle properties
        E, px, py, pz, ityp, chg = event.T
        p = np.sqrt(px**2.+py**2.+pz**2.)
        pt = np.sqrt(px**2.+py**2.)
        eta = 0.5*np.log((p+pz)/(p-pz))

        # cut on rapidity and charge
        cut_pi = (np.abs(eta) < 1.0) & (ityp.astype(int) == 101) & (chg.astype(int) != 0)
        cut_k  = (np.abs(eta) < 1.0) & (ityp.astype(int) == 106) & (chg.astype(int) != 0)
        cut_p  = (np.abs(eta) < 1.0) & (ityp.astype(int) == 1) & (chg.astype(int) == 1)

        # append pt
        pion.extend(pt[cut_pi])
        kaon.extend(pt[cut_k])
        proton.extend(pt[cut_p])

        # increment event counter
        nev += 1
        print(nev)

    pion = np.array(pion)
    kaon = np.array(kaon)
    proton = np.array(proton)
    arrays = pion, kaon, proton
    names = "/pion", "/kaon", "/proton"

    for array, name in zip(arrays, names):
        
        hist, bin_edges = np.histogram(array, bins=np.linspace(0., 3., 20), weights=1/(nev*2*np.pi*array), normed=False)
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

        directory="results/" + eos + centrality
        if not os.path.exists(directory):
            os.makedirs(directory)

        np.savez(directory + name, hist=hist, bin_centers=bin_centers)
