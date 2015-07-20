#!/usr/bin/env python3


from __future__ import division, print_function, unicode_literals

import functools

import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib import cm
import numpy as np
import seaborn as sns
import os, sys

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
    centralities = '/0to5/', '/20to30/', '/40to50/'
    labels = '0–5%', '20–30%', '40–50%'
    data_columns = 3, 11, 15
    
    fig = plt.figure(figsize=(textwidth, textwidth*aspect))
    gs = plt.GridSpec(3, 9, wspace=0.0, hspace=0.0)
    ax1 = fig.add_subplot(gs[0:1,0:3])
    ax2 = fig.add_subplot(gs[0:1,3:6], sharey=ax1)
    ax3 = fig.add_subplot(gs[0:1,6:9], sharey=ax1)
    axes = ax1, ax2, ax3

    for ax in axes:
        despine(ax)
        ax.tick_params(top='off', right='off', labelbottom='off', labelleft='off')
    ax1.tick_params(labelleft='on')

    for ax, centrality, label, data_column in zip(axes, centralities, labels, data_columns):
        try:
            # load model results
            pi = np.load("../results/HotQCD/" + centrality + "/pion.npz")
            k = np.load("../results/HotQCD/" + centrality + "/kaon.npz")
            p = np.load("../results/HotQCD/" + centrality + "/proton.npz")
            
            # load exp results
            pi_pt, pi_yield, pi_yield_err = np.loadtxt("exp/pi+.dat", usecols=(0, data_column, data_column + 1), unpack=True)
            k_pt, k_yield, k_yield_err = np.loadtxt("exp/k+.dat", usecols=(0, data_column, data_column + 1), unpack=True)
            p_pt, p_yield, p_yield_err = np.loadtxt("exp/p.dat", usecols=(0, data_column, data_column + 1), unpack=True)

        except IOError:
            print("missing file in", system + centrality)
        else:
            # plot model results
            ax.plot(pi['bin_centers'], 100*pi['hist'], linewidth=0.75, color=plt.cm.Blues(0.6), label=r'$\pi^+ \! \times 10^2$')
            ax.fill_between(pi['bin_centers'], 100*(pi['hist'] - 2*pi['err']), 100*(pi['hist'] + 2*pi['err']), color=plt.cm.Blues(0.6), alpha=0.4, linewidth=0)

            ax.plot(k['bin_centers'], 10*k['hist'], linewidth=0.75, color=plt.cm.Oranges(0.6), label=r'$K^+ \! \times 10$')
            ax.fill_between(k['bin_centers'], 10*(k['hist'] - 2*k['err']), 10*(k['hist'] + 2*k['err']), color=plt.cm.Oranges(0.6), alpha=0.4, linewidth=0)

            ax.plot(p['bin_centers'], p['hist'], linewidth=0.75, color=plt.cm.Greens(0.6), label=r'$p$')
            ax.fill_between(p['bin_centers'], (p['hist'] - 2*p['err']), (p['hist'] + 2*p['err']), color=plt.cm.Greens(0.6), alpha=0.4, linewidth=0)

            # plot exp results
            ax.plot(pi_pt, 100*pi_yield, 'o', color=offblack, markerfacecolor=offblack, markeredgecolor='none', markersize=3, zorder=100, label=r'$\pi^+ \! \times 10^2$')
            ax.plot(k_pt, 10*k_yield, 's', color=offblack, markerfacecolor=offblack, markeredgecolor='none', markersize=2.6, zorder=100, label=r'$K^+ \! \times 10$')
            ax.plot(p_pt, p_yield, '^', color=offblack, markerfacecolor=offblack, markeredgecolor='none', markersize=3.5, zorder=100, label=r'$p$')

            # x-axis labels
            ax.set_xlim([0,3.0])
            ax.set_xticks([0.0, 0.5,  1.0, 1.5,  2.0, 2.5])
        
            # y-axis labels
            ax.set_ylim([10e-4, 10e4])
            ax.set_yscale('log')
            ax.minorticks_off()
        
        # annotations
        ax.annotate(label, xy=(0.5,0.95), xycoords='axes fraction', ha='center', fontsize=texnormal)
   
    handles, labels = ax1.get_legend_handles_labels()
    ax2.legend(handles[0:3], labels[0:3], bbox_to_anchor=(1.05,1.05), handletextpad=0.5, handlelength=1.4, labelspacing=0)
    ax3.annotate('PHENIX', xy=(0.95,1), xycoords='axes fraction', ha='center')
    ax3.legend(handles[3:6], labels[3:6], bbox_to_anchor=(1.15,1.05), handletextpad=0.0, labelspacing=0)

    # y-axis labels
    ax1.set_ylabel(r'$dN/(2 \pi p_T dp_T dy) \,[\mathrm{GeV}^{-2}]$')
    ax3.annotate('HotQCD', xy=(1.025,0.5), va='center', xycoords='axes fraction', rotation=-90, fontsize=texnormal)
    ax1.set_yticks([10e-4, 10e-2, 10e0, 10e2])
    
    # WB/HotQCD
    ax1a = fig.add_subplot(gs[1,0:3], sharex=ax1)
    ax2a = fig.add_subplot(gs[1,3:6], sharey=ax1a)
    ax3a = fig.add_subplot(gs[1,6:9], sharey=ax1a)
    axes = ax1a, ax2a, ax3a

    for ax in axes:
        despine(ax)
        ax.tick_params(top='off', right='off', labelbottom='off', labelleft='off')
    ax1a.tick_params(labelleft='on')
    ax1a.set_ylabel("Ratio")
    ax1a.set_yticks([0.7, 1.0, 1.3])
    ax3a.annotate('WB', xy=(1.025,0.5), va='center', xycoords='axes fraction', rotation=-90, fontsize=texnormal)

    for ax, centrality in zip(axes, centralities):
        try:
            pi_HotQCD = np.load("../results/HotQCD/" + centrality + "/pion.npz")
            k_HotQCD = np.load("../results/HotQCD/" + centrality + "/kaon.npz")
            p_HotQCD = np.load("../results/HotQCD/" + centrality + "/proton.npz")
            pi_WB = np.load("../results/WB/" + centrality + "/pion.npz")
            k_WB = np.load("../results/WB/" + centrality + "/kaon.npz")
            p_WB = np.load("../results/WB/" + centrality + "/proton.npz")
        except IOError:
            print("missing file in", system + centrality)
        else:
            pt = pi_HotQCD['bin_centers']
            ratio = pi_WB['hist']/pi_HotQCD['hist']
            err = ratio*np.sqrt((2*pi_WB['err']/pi_WB['hist'])**2 + (2*pi_HotQCD['err']/pi_HotQCD['hist'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Blues(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Blues(0.6), alpha=0.4, linewidth=0)  

            pt = k_HotQCD['bin_centers']
            ratio = k_WB['hist']/k_HotQCD['hist']
            err = ratio*np.sqrt((2*k_WB['err']/k_WB['hist'])**2 + (2*k_HotQCD['err']/k_HotQCD['hist'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Oranges(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Oranges(0.6), alpha=0.4, linewidth=0)

            pt = p_HotQCD['bin_centers']
            ratio = p_WB['hist']/p_HotQCD['hist']
            err = ratio*np.sqrt((2*p_WB['err']/p_WB['hist'])**2 + (2*p_HotQCD['err']/p_HotQCD['hist'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Greens(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Greens(0.6), alpha=0.4, linewidth=0)

            ax.plot(np.linspace(0,3,100), np.ones(100), linewidth=0.2, color='gray')
           
            ax.set_xlim([0,3.0])
            ax.set_ylim([0.5,1.5])

    ax2a.tick_params(labelleft='off')
    ax3a.tick_params(labelleft='off')

    # s95-PCE/HotQCD
    ax1b = fig.add_subplot(gs[2,0:3], sharex=ax1)
    ax2b = fig.add_subplot(gs[2,3:6], sharex=ax1, sharey=ax1b)
    ax3b = fig.add_subplot(gs[2,6:9], sharex=ax1, sharey=ax1b)
    axes = ax1b, ax2b, ax3b

    for ax in axes:
        despine(ax)
        ax.tick_params(top='off', right='off', labelleft='off')
    ax1b.tick_params(labelleft='on')
    ax1b.set_ylabel("Ratio")
    ax1b.set_yticks([0.7, 1.0, 1.3])
    ax3b.annotate('S95', xy=(1.025,0.5), va='center', xycoords='axes fraction', rotation=-90, fontsize=texnormal)

    for ax, centrality in zip(axes, centralities):
        try:
            pi_HotQCD = np.load("../results/HotQCD/" + centrality + "/pion.npz")
            k_HotQCD = np.load("../results/HotQCD/" + centrality + "/kaon.npz")
            p_HotQCD = np.load("../results/HotQCD/" + centrality + "/proton.npz")
            pi_s95 = np.load("../results/s95/" + centrality + "/pion.npz")
            k_s95 = np.load("../results/s95/" + centrality + "/kaon.npz")
            p_s95 = np.load("../results/s95/" + centrality + "/proton.npz")
        except IOError:
            print("missing file in", system + centrality)
        else:

            pt = pi_HotQCD['bin_centers']
            ratio = pi_s95['hist']/pi_HotQCD['hist']
            err = ratio*np.sqrt((2*pi_s95['err']/pi_s95['hist'])**2 + (2*pi_HotQCD['err']/pi_HotQCD['hist'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Blues(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Blues(0.6), alpha=0.4, linewidth=0)  

            pt = k_HotQCD['bin_centers']
            ratio = k_s95['hist']/k_HotQCD['hist']
            err = ratio*np.sqrt((2*k_s95['err']/k_s95['hist'])**2 + (2*k_HotQCD['err']/k_HotQCD['hist'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Oranges(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Oranges(0.6), alpha=0.4, linewidth=0)

            pt = p_HotQCD['bin_centers']
            ratio = p_s95['hist']/p_HotQCD['hist']
            err = ratio*np.sqrt((2*p_s95['err']/p_s95['hist'])**2 + (2*p_HotQCD['err']/p_HotQCD['hist'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Greens(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Greens(0.6), alpha=0.4, linewidth=0)

            ax.plot(np.linspace(0,3,100), np.ones(100), linewidth=0.2, color='gray')
            
            ax.set_xlim([0,3.0])
            ax.set_xlabel("$p_T$ [GeV]")
            ax.set_ylim([0.5,1.5])

    plt.tight_layout(pad=.1, w_pad=0, rect=[0,0,0.975,1])

if __name__== "__main__":
    spectra()
