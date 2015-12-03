#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
import random 

charged = np.sort([np.load(eos)['charged'] for eos in sys.argv[1:]])
charged_err = np.sort([np.load(eos)['charged_err'] for eos in sys.argv[1:]])
plt.errorbar((np.ones(len(charged)) + 0.04*np.random.rand(len(charged))), charged, \
            yerr=charged_err, fmt='o', color=plt.cm.Blues(0.65), capsize=3, lw=2, ms=4)
plt.xlim(0,2)
#plt.ylim(0.56, 0.58)
plt.show()
