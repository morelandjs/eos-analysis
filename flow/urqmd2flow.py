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
            event.append((E, px, py, pz, ityp)) 
        else:
            if event:
                yield np.array(event, dtype=float)
                event.clear()

    if event: 
        yield np.array(event, dtype=float)
                            

for event in events():
    E, px, py, pz, ityp = event.T
    phi = np.arctan2(py, px)
    pt = np.sqrt(px**2.+py**2.) 

    for i in zip(phi, pt, ityp):
        print(*i)
    print()
