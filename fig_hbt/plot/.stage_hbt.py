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


# construct the Gaussian fit function
#def residuals(args, A): 
#    norm, lmbda, Ro, Rs, Rl = args
#    qo, qs, ql = np.mgrid[-0.5:0.5:20j,-0.5:0.5:20j,-0.5:0.5:20j]
#    B = norm*(1. + lmbda * np.exp(-(qo*Ro)**2. - (qs*Rs)**2. - (ql*Rl)**2.))
#    return np.array(A - B).flatten()

# Gaussian fit function
def gaussian(q, norm, lmbda, Ro, Rs, Rl):
    return np.ravel(norm*(1. + lmbda * np.exp(-(q[0]*Ro)**2. - (q[1]*Rs)**2. - (q[2]*Rl)**2.)))

# set path
path = sys.argv[1]

# specify EOS folder, e.g. path/HotQCD
#for centrality in "/0to10/", "/10to20/", "/20to40/":
for centrality in "/0to10/",:

    # pull eos name from end of path
    eos = os.path.basename(sys.argv[1])
    print(centrality)

    # initialize correlation containers
    Csame = np.zeros((10,20,20,20))
    Cmixed = np.zeros((10,20,20,20))

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

    # loop over kt bins
    kt_list = np.linspace(0,1,10)
    results = []
    errors = []
    for ikt, kt in enumerate(kt_list):
        q = np.array(np.mgrid[-0.5:0.5:20j,-0.5:0.5:20j,-0.5:0.5:20j])
        popt, pcov = scipy.optimize.curve_fit(gaussian, q, np.ravel(C[ikt]), p0 = [1.0, 0.4, 4, 4, 4]) 
        err = np.diagonal(np.sqrt(pcov))
        
        #leastsq = scipy.optimize.leastsq(residuals, p0, args=(C[ikt]), full_output=1)
        #vals = leastsq[0]
        #cov = leastsq[1]*np.var(residuals(vals, C[ikt]))
        #err = np.sqrt(np.array((cov[0,0], cov[1,1], cov[2,2], cov[3,3], cov[4,4])))
        
        results.append(popt)
        errors.append(err)

    #norm, lmbda, Ro, Rs, Rl = np.array(np.abs(results)).T
    #hbt_radii = Ro, Rs, Rl

    #plt.plot(q,  norm[0]*(1. + lmbda[0] * np.exp(-(q*Rl[0])**2.)))
    #plt.show()
    #exit()

#    directory="../results-staged/" + eos + centrality
#    if not os.path.exists(directory):
#        os.makedirs(directory)

#    np.savez(directory + "R", kt_list=kt_list, Ro=Ro, Rs=Rs, Rl=Rl)
