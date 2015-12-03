#!/usr/bin/env python


from __future__ import division, print_function, unicode_literals

import functools
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib import cm
import numpy as np
import glob, os, sys
import h5py
import scipy
from scipy.optimize import leastsq


# Gaussian fit function
def gaussian(q, norm, lmbda, Ro, Rs, Rl):
    hbarc = 0.19733
    return np.ravel(norm*(1. + lmbda * np.exp(-(q[0]*Ro/hbarc)**2 -(q[1]*Rs/hbarc)**2 -(q[2]*Rl/hbarc)**2)))

# specify EOS folder, e.g. path/HotQCD
for centrality in "/0to10/", "/10to20/", "/20to40/":
#for centrality in "/0to10/",:

    # pull eos name from end of path
    path = sys.argv[1]
    eos = os.path.basename(path)
    print(centrality)

    # initialize correlation containers
    Csame = np.zeros((10,15,15,15))
    Cmixed = np.zeros((10,15,15,15))

    # calculate same event correlations
    for fsame, fmixed in zip(os.listdir(path + centrality + "/same"), os.listdir(path + centrality + "/mixed")):
        
        # open same event
        file_path = path + centrality + "/same/" + fsame
        file = h5py.File(file_path, 'r')
        Csame += np.nan_to_num(np.array(file['correlations'][()]))/float(len(os.listdir(path + centrality + "/same")))
        file.close()
        
        # open mixed event
        file_path = path + centrality + "/mixed/" + fmixed
        file = h5py.File(file_path, 'r', fillValue=0.)
        Cmixed += np.nan_to_num(np.array(file['correlations'][()]))/float(len(os.listdir(path + centrality + "/mixed")))
        file.close()

    # calculate the correlation matrix C = Csame/Cmixed
    C = np.true_divide(Csame, Cmixed)
    #plt.imshow(C[6,:,9,:])
    #plt.show()
    #exit()

    # initialize containers
    kt_list = np.linspace(0,1,10)
    results = []
    errors = []
    
    # loop over kt bins
    for ikt, kt in enumerate(kt_list):
        q = np.array(np.mgrid[0.0:0.1:15j,0.0:0.1:15j,0.0:0.1:15j])
        popt, pcov = scipy.optimize.curve_fit(gaussian, q, np.ravel(C[ikt])) 
        results.append(np.abs(popt))
        errors.append(np.diagonal(np.sqrt(np.abs(pcov))))

    norm, lmbda, Ro, Rs, Rl = np.array(results).T
    dnorm, dlmbda, dRo, dRs, dRl = np.array(errors).T

    # write to file
    directory="../results-staged/" + eos + centrality
    if not os.path.exists(directory):
        os.makedirs(directory)
    np.savez(directory + "R", kt_list=kt_list, Ro=Ro, dRo=dRo, Rs=Rs, dRs=dRs, Rl=Rl, dRl=dRl)

    # plot correlation function vs fit
    #resolution = 72.27
    #textwidth = 510/resolution
    #fig = plt.figure(figsize=(textwidth, textwidth/4))
    #gs = plt.GridSpec(1, 3, wspace=0.0, hspace=0.0)
    
    #q = np.linspace(0.0,0.1,15)
    #kt_indices = 1, 3, 5, 7, 9
    #kt_indices = 4,
    #Cslices = C[:,:,9,9
    #Rslices = Ro, Rs, Rl
    #for col, Cslice, Rslice in zip(range(3), Cslices, Rslices):
    #    ax = fig.add_subplot(gs[0,col])
    #    for kt_index in kt_indices:
    #        ax.plot(q, Cslice[kt_index,:], 'o')
    #        hbarc = 0.19733
    #        ax.plot(q, norm[kt_index]*(1 + lmbda[kt_index]*np.exp(-(q*Rslice[kt_index]/hbarc)**2)), 'r--')
    #plt.show()
