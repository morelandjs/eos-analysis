#!/usr/bin/env python3


from __future__ import division, print_function, unicode_literals

import functools
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib import cm
import numpy as np
import os, sys, glob

# figure properties
aspect = 1/1.618
resolution = 72.27
columnwidth = 246/resolution
textwidth = 510/resolution
textiny, texsmall, texnormal = 8., 9.25, 10.
offblack = '#262626'

plt.rcdefaults()
plt.rcParams.update({
        'font.family': ['CMU Serif'],
        'font.size': texsmall,
        'legend.fontsize': 8.5,
        'axes.labelsize': texsmall,
        'axes.titlesize': texsmall,
        'xtick.labelsize': textiny,
        'ytick.labelsize': textiny,
        'lines.linewidth': .9,
        'patch.linewidth': .9,
        'axes.linewidth': .5,
        'text.color': offblack,
        'axes.edgecolor': offblack,
        'axes.labelcolor': offblack,
        'xtick.color': offblack,
        'ytick.color': offblack,
        'xtick.major.size': 2,
        'ytick.major.size': 2,
        'xtick.major.pad': 1.8,
        'ytick.major.pad': 1.8,
        'legend.numpoints': 1,
        'legend.scatterpoints': 1,
        'legend.frameon': False,
        'pdf.fonttype': 42
})

def plotfn(f):
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        print(f.__name__)
        f(*args, **kwargs)
        plt.savefig('{}.pdf'.format(f.__name__))
        plt.close()
    return wrapper


def despine(ax=None, remove_ticks=False):
    if ax is None:
        ax = plt.gca()
    for spine in 'top', 'right':
        ax.spines[spine].set_visible(False)
    if remove_ticks:
        for ax_name in 'xaxis', 'yaxis':
            getattr(ax, ax_name).set_ticks_position('none')
    else:
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')


@plotfn
def spectra():
    fig = plt.figure(figsize=(columnwidth, 1.2*columnwidth))
    gs = plt.GridSpec(5, 1, wspace=0.0, hspace=0.0)
    ax1 = fig.add_subplot(gs[0:2,0])
    ax2 = fig.add_subplot(gs[2,0], sharex=ax1)
    ax3 = fig.add_subplot(gs[3,0], sharex=ax1)
    ax4 = fig.add_subplot(gs[4,0], sharex=ax1)
    axes = ax1, ax2, ax3, ax4

    for ax in axes:
        despine(ax)
        ax.tick_params(top='off', right='off', labelbottom='off', labelleft='on')
    ax1.tick_params(top='off', right='off', labelbottom='off', labelleft='on')
    ax4.tick_params(top='off', right='off', labelbottom='on', labelleft='on')

    splines = glob.glob('../results/*')
    hqcd = '/home/jsm55/Research/eos-analysis/spectra/results/hotqcd/'

    particles = 'pion', 'kaon', 'proton'
    labels = r'$\pi^+ \times 10^2$', r'$K^+ \! \times 10$', r'$p$'
    offsets = 100, 10, 1
    colors = plt.cm.Blues(0.6), plt.cm.Oranges(0.6), plt.cm.Greens(0.6)
    paths = "exp/pi+.dat", "exp/k+.dat", "exp/p.dat"
    markers = 'o', 's', '^'
    sizes = 3, 2.6, 3.5
    rlabels = r'$\pi^+ \! / \langle \pi^+ \! \rangle$', r'$K^+ \! / \langle K^+ \! \rangle$', r'$p/\langle p \rangle$' 

    # loop over particle types
    for particle, offset, color, ax, label, path, marker, size, rlabel in zip(particles, offsets, colors, axes[1:4], labels, paths, markers, sizes, rlabels):

        # Experimental yield plots
        x, y , dy = np.loadtxt(path, usecols=(0, 11, 12), unpack=True)
        ax1.errorbar(x, y*offset, yerr=dy*offset, fmt=marker, capsize=0, markersize=size, color=offblack, label=label)

        # Yield plots
        pt = np.array([np.load(spline + "/20to30/" + particle + ".npz")['bin_centers'] for spline in splines])
        hist_sample = np.array([np.load(spline + "/20to30/" + particle + ".npz")['hist'] for spline in splines])
        err_sample = np.array([np.load(spline + "/20to30/" + particle + ".npz")['err'] for spline in splines])
        hist_hqcd = np.array(np.load(hqcd + "/20to30/" + particle + ".npz")['hist'])
        err_hqcd = np.array(np.load(hqcd + "/20to30/" + particle + ".npz")['err'])
        [ax1.plot(pt[index], offset*hist_sample[index], linewidth=0.2, color=color, label=label if index==0 else "") for index in range(len(splines))]
        
        # Ratio plots
        ratio = [hist_sample[index]/hist_hqcd for index in range(len(splines))]
        stat_err = [ratio[index] * np.sqrt((2*err_sample[index]/hist_sample[index])**2 + (2*err_hqcd/hist_hqcd)**2) for index in range(len(splines))]
        sys_err = np.percentile(ratio, [2.5, 97.5], axis=0) 
        ax.plot(np.mean(pt, axis=0), sys_err[0], dashes=(4,1), linewidth=0.6, color=color)
        ax.plot(np.mean(pt, axis=0), sys_err[1], dashes=(4,1), linewidth=0.6, color=color)
        ax.fill_between(np.mean(pt, axis=0), np.mean(ratio, axis=0) - np.mean(stat_err, axis=0), np.mean(ratio, axis=0) + np.mean(stat_err, axis=0), color=color, alpha=0.4, linewidth=0)
        ax.annotate(rlabel, xy=(0.97,0.5), va='center', xycoords='axes fraction', rotation=-90)
         
    # ax1 properties        
    ax1.set_xlim([0,3.0])
        
    ax1.set_xticks([0.0, 0.5,  1.0, 1.5,  2.0, 2.5])
    ax1.set_yscale('log')
    ax1.set_ylim([10e-4, 10e4])
    ax1.set_yticks([1e-3, 1e-1, 1e1, 1e3])
    ax1.set_ylabel(r'$dN/(2 \pi p_T dp_T dy) \,[\mathrm{GeV}^{-2}]$')
    ax1.minorticks_off()
    ax1.annotate('HQCD samples', xy=(0.775,1), xycoords='axes fraction', va='top', ha='right')
    ax1.annotate('PHENIX', xy=(1.03,1), xycoords='axes fraction', va='top', ha='right')
    ax1.annotate('20–30% centrality', xy=(0.05,0.08), xycoords='axes fraction', va='bottom', ha='left')
    
    # ax2-ax4 properties
    ax2.set_ylim([0.925, 1.075])
    ax2.set_yticks([0.95, 1, 1.05])
    ax2.set_yticklabels([0.95, 1, 1.05])
    ax2.set_ylabel('Ratio')
    
    ax3.set_ylim([0.925, 1.075])
    ax3.set_yticks([0.95, 1, 1.05])
    ax3.set_yticklabels([0.95, 1, 1.05])
    ax3.set_ylabel('Ratio')

    ax4.set_ylim([0.925, 1.075])
    ax4.set_yticks([0.95, 1, 1.05])
    ax4.set_yticklabels([0.95, 1, 1.05])
    ax4.set_ylabel('Ratio')
    ax4.set_xlabel(r'$p_T$ [GeV]')

    # legend
    handles, labels = ax1.get_legend_handles_labels()
    leg=ax1.legend(handles, labels, bbox_to_anchor=(1.09, 0.96), ncol=2, columnspacing=0.1, handletextpad=0.5, handlelength=1.4, labelspacing=0)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(1.0)    
    plt.tight_layout(pad=0.1, w_pad=0, rect=[0,0,0.975,1])

if __name__== "__main__":
    spectra()
