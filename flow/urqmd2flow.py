#!/usr/bin/env python
import numpy as np
import math
import sys
import gzip

# parse spectra (avg over events)
header=True
E,px,py,pz,ityp = 0.,0.,0.,0.,0
for fn in sys.argv[1:]: 
    for l in gzip.open(fn):
        try:
	    p = l[97:192].replace(b'D',b'E')
	    E = float(p[0:23])
            px = float(p[24:47])
            py = float(p[48:71])
	    pz = float(p[72:95])
            ityp = int(l[218:221])
	except ValueError:
            if not header:
                header = True
		sys.stdout.write("\n")
	else:
	    if header:	
		header=False

    	    phi = np.arctan2(py,px)
    	    pt = np.sqrt(px**2.+py**2.)
    	    y = 0.5*np.log((E+pz)/max(E-pz,1e-12))
	    row = (ityp, pt, phi)

            # print output to file
            sys.stdout.write("{:20}{:20}{:20}\n".format(*row))
