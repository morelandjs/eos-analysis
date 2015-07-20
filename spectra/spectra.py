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
#for centrality in "/0to5/", "/20to30/", "/40to50/":
for centrality in "/20to30/", "/40to50/":

    # max eta cut
    eta_max = 1.0

    eos = os.path.basename(sys.argv[1])
    filenames = glob.glob(sys.argv[1] + centrality + "/*.urqmd.gz")
    print(centrality)

    pion = []
    kaon = []
    proton = []
    nch = []
    nev = 0

    # loop over events in each folder
    for event in itertools.chain.from_iterable(events(f) for f in filenames[:]):

        # unpack particle properties
        E, px, py, pz, ityp, chg = event.T
        p = np.sqrt(px**2.+py**2.+pz**2.)
        pt = np.sqrt(px**2.+py**2.)
        rapidity = 0.5*np.log((p+pz)/(p-pz))
        #rapidity = 0.5*np.log((E+pz)/(E-pz))

        # cut on rapidity and charge
        cut_pi  = (np.abs(rapidity) < eta_max) & (ityp.astype(int) == 101) & (chg.astype(int) > 0)
        cut_k   = (np.abs(rapidity) < eta_max) & (ityp.astype(int) == 106) & (chg.astype(int) > 0)
        cut_p   = (np.abs(rapidity) < eta_max) & (ityp.astype(int) == 1) & (chg.astype(int) == 1)
        cut_nch = (np.abs(rapidity) < eta_max) & (chg.astype(int) != 0)

        # append pt
        pion.extend(pt[cut_pi])
        kaon.extend(pt[cut_k])
        proton.extend(pt[cut_p])
        nch.append(pt[cut_nch].size)

        # increment event counter
        nev += 1

    pion = np.array(pion)
    kaon = np.array(kaon)
    proton = np.array(proton)
    nch = np.array(nch)

    # print particle density
    print("charged:", nch.mean()/(2.*eta_max))
    print("pion:", pion.size/float(nev)/(2.*eta_max))
    print("kaon:", kaon.size/float(nev)/(2.*eta_max))
    print("proton:", proton.size/float(nev)/(2.*eta_max))
    
#    pt_lists = pion, kaon, proton
#    names = "/pion", "/kaon", "/proton"

#    for pt_list, name in zip(pt_lists, names):
#        pt_min = 0. 
#        pt_max = 3.2
#        dpt = 0.2 
#        deta = 2.*eta_max
#        bin_edges = np.arange(pt_min, pt_max, dpt)
#        bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
#        w = 2.*np.pi*pt_list*dpt*deta*nev
#        hist = np.histogram(pt_list, bins=bin_edges, weights=1./w, normed=False)[0]
#        err = np.sqrt(np.histogram(pt_list, bins=bin_edges, normed=False)[0])/(2.*np.pi*bin_centers*dpt*deta*nev)

#        directory="results/" + eos + centrality
#        if not os.path.exists(directory):
#            os.makedirs(directory)
#
#        np.savez(directory + name, bin_centers=bin_centers, hist=hist, err=err)
