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
            if len(l) == 308:
                p = l[65:128].replace(b'D', b'E')
                E = p[0:15]
                px = p[16:31]
                py = p[32:47]
                pz = p[48:63]
                ityp = l[145:155]
                chg = l[159:161]
                event.append((E, px, py, pz, ityp, chg))
            else:
                if event:
                    try:
                        yield np.array(event, dtype=float)
                    except ValueError:
                        print('filename: ', filename, 'line: ', l)
                    event.clear()

        if event:
            yield np.array(event, dtype=float)


# max eta cut
eta_max = 0.5
centrality = "20to30"

# loop over equations of state
for eos, path in zip([os.path.basename(path) for path in sys.argv[1:]], sys.argv[1:]):

    # progress indicator
    print(eos)

    # initialize arrays
    pion = []
    kaon = []
    proton = []
    charged = []

    # loop over events in each eos folder
    for event in itertools.chain.from_iterable(events(f) for f in glob.glob(path + "/" + centrality + "/*.urqmd.gz")):

        # unpack particle properties
        E, px, py, pz, ityp, chg = event.T
        p = np.sqrt(px**2.+py**2.+pz**2.)
        pt = np.sqrt(px**2.+py**2.)
        rapidity = 0.5*np.log((p+pz)/(p-pz))

        # cut on rapidity and charge
        cut_pi  = (np.abs(rapidity) < eta_max) & (ityp.astype(int) == 101) & (chg.astype(int) > 0)
        cut_k   = (np.abs(rapidity) < eta_max) & (ityp.astype(int) == 106) & (chg.astype(int) > 0)
        cut_p   = (np.abs(rapidity) < eta_max) & (ityp.astype(int) == 1) & (chg.astype(int) == 1)
        cut_chg = (np.abs(rapidity) < eta_max) & (chg.astype(int) != 0)

        # append the mean pt
        if pt[cut_pi].size > 0: pion.append(pt[cut_pi].mean()) 
        if pt[cut_k].size > 0: kaon.append(pt[cut_k].mean()) 
        if pt[cut_p].size > 0: proton.append(pt[cut_p].mean()) 
        if pt[cut_chg].size > 0: charged.append(pt[cut_chg].mean()) 

    # convert to array format
    pion = np.array(pion)
    kaon = np.array(kaon)
    proton = np.array(proton)
    charged = np.array(charged)

    # save to file
    directory = "results/samples/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    np.savez(directory + eos + ".npz", \
            pion=pion.mean(), pion_err=pion.std()/np.sqrt(pion.size), \
            kaon=kaon.mean(), kaon_err=kaon.std()/np.sqrt(kaon.size), \
            proton=proton.mean(), proton_err=proton.std()/np.sqrt(proton.size), \
            charged=charged.mean(), charged_err=charged.std()/np.sqrt(charged.size) \
            )

