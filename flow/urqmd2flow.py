#!/usr/bin/env python3

import numpy as np
import sys


def events():
    event = []
    for l in sys.stdin.buffer:
        if len(l) > 200: 
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
                            

for event in events():

    # unpack particle properties
    E, px, py, pz, ityp, chg = event.T
    phi = np.arctan2(py, px)
    p = np.sqrt(px**2.+py**2.+pz**2.)
    pt = np.sqrt(px**2.+py**2.)
    eta = 0.5*np.log((p+pz)/(p-pz))

    # cut on rapidity and charge
    cut = (np.abs(eta) < 0.5) & (chg.astype(int) != 0)

    for i in zip(ityp[cut], pt[cut], phi[cut]):
        print(*i)
    print()
