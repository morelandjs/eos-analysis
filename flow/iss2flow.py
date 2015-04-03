#!/usr/bin/env python3

import numpy as np
import sys

pt, v2_, v3_ = np.array([np.loadtxt(f, usecols=(0,8,11)) for f in sys.argv[1:]]).transpose(2,0,1)
pt = pt[0]
v2 = np.mean(v2_, axis=0)
dv2 = np.std(v2_, axis=0)/np.sqrt(v2_.size)
v3 = np.mean(v3_, axis=0)
dv3 = np.std(v3_, axis=0)/np.sqrt(v3_.size)

for pt, v2, dv2, v3, dv3 in zip(pt, v2, dv2, v3, dv3):
   print(pt, v2, dv2, v3, dv3)
