#!/usr/bin/env python
import os
import sys
import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy
import seaborn as sns
from scipy.optimize import leastsq

# seaborn properties
sns.set(style="white", font='CMU Serif', font_scale=1.5, rc={"lines.linewidth": 2.5})

# initialize correlation containers
Csame = np.zeros((10,20,20,20))
Cmixed = np.zeros((10,20,20,20))

# construct the Gaussian fit function
def residuals(args, A): 
    lmbda, Ro, Rs, Rl = args
    qo, qs, ql = np.mgrid[-0.5:0.5:20j,-0.5:0.5:20j,-0.5:0.5:20j]
    B = 1. + lmbda * np.exp(-(qo*Ro)**2. - (qs*Rs)**2. - (ql*Rl)**2.)
    err = np.array(A - B).flatten()
    return err

# create the figure template
fig, axes = plt.subplots(ncols=3, figsize=(4.6, 3), sharey=True)

# folder containing hbt data
for folder in sys.argv[1:]:
    
    # print working folder
    print(folder)

    # calculate same event correlations
    for filename in os.listdir(folder + "/same"):
        file_path = folder + "/same/" + filename
        file = h5py.File(file_path, 'r')
        Csame += np.asarray(file['correlations'][()])/float(len(os.listdir(folder)))
        file.close()

    # calculate mixed event correlations
    for filename in os.listdir(folder + "/mixed"):
        file_path = folder + "/mixed/" + filename
        file = h5py.File(file_path, 'r')
        Cmixed += np.asarray(file['correlations'][()])/float(len(os.listdir(folder)))
        file.close()

    # calculate the correlation matrix C = Csame/Cmixed
    C = np.true_divide(Csame, Cmixed)

    # loop over kt bins
    results = []
    kt_list = np.linspace(0,1,10)
    for ikt, kt in enumerate(kt_list):
        A = C[ikt]
        p0 = [0.5, 6, 6, 9]
        results.append(scipy.optimize.leastsq(residuals, p0, args=(A))[0])

    lmbda, Ro, Rs, Rl = np.array(np.abs(results)).T
    hbt_radii = Ro, Rs, Rl
    labels = "$R_\t{out}$", "$R_\t{side}$", "$R_\t{long}$"
    for ax, R, label in zip(axes, hbt_radii, labels):
        ax.plot(kt_list, R, linewidth=0.5)
        ax.set_ylim([2,11.1])
        ax.set_xticks(np.arange(0.0,1,0.5))
        ax.set_yticks(np.arange(2, 11, 2))
        ax.annotate(label, xy=(0.35,10), xycoords='data', fontsize=16)
        ax.set_xlabel("kt [GeV]")

axes[0].set_ylabel("R [fm]")
plt.tight_layout(pad=0, w_pad=0)
plt.show()
