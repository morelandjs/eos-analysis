#!/usr/bin/env python3


from __future__ import division, print_function, unicode_literals

import functools
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib import cm
import numpy as np
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
def v2():
    centralities = '/0to10/', '/20to30/', '/40to50/'
    labels = '0–10%', '20–30%', '40–50%'
    
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

    for ax, centrality, label in zip(axes, centralities, labels):
        try:
            pi = np.load("../results/hotqcd/" + centrality + "/pion.npz")
            k = np.load("../results/hotqcd/" + centrality + "/kaon.npz")
            p = np.load("../results/hotqcd/" + centrality + "/proton.npz")
        except IOError:
            print("missing file in", system + centrality)
        else:
            ax.plot(pi['pt'], pi['v2'], linewidth=0.75, color=plt.cm.Blues(0.6), label='$\pi^\pm$')
            ax.fill_between(pi['pt'], pi['v2'] - pi['dv2'], pi['v2'] + pi['dv2'], color=plt.cm.Blues(0.6), alpha=0.5, linewidth=0)
           
            ax.plot(k['pt'], k['v2'], linewidth=0.75, color=plt.cm.Oranges(0.6), label='$K^\pm$')
            ax.fill_between(k['pt'], k['v2'] - k['dv2'], k['v2'] + k['dv2'], color=plt.cm.Oranges(0.6), alpha=0.5, linewidth=0)

            ax.plot(p['pt'], p['v2'], linewidth=0.75, color=plt.cm.Greens(0.6), label='$p$')
            ax.fill_between(p['pt'], p['v2'] - p['dv2'], p['v2'] + p['dv2'], color=plt.cm.Greens(0.6), alpha=0.5, linewidth=0)

        # x-axis labels
        ax.set_xlim([0,2.5])
        ax.set_xticks([0.0, 0.5,  1.0, 1.5,  2.0])
        
        # y-axis labels
        ax.set_ylim([0., 0.175])
        ax.minorticks_off()
        
        # annotations
        ax.annotate(label, xy=(0.5,0.9), xycoords='axes fraction', ha='center', fontsize=texnormal)
    
    ax3.legend(bbox_to_anchor=(0.35,1.05), handletextpad=0.5, handlelength=2.0, labelspacing=0)

    # y-axis labels
    ax1.set_ylabel(r'$v_2$')
    ax3.annotate('HQCD', xy=(1.025,0.5), va='center', xycoords='axes fraction', rotation=-90, fontsize=texnormal)
    ax1.set_yticks([0., 0.05, 0.1, 0.15])
    
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
    ax1a.set_yticks([0.9, 1.0, 1.1])
    ax3a.annotate('WB/HQCD', xy=(1.025,0.5), va='center', xycoords='axes fraction', rotation=-90, fontsize=texnormal)

    for ax, centrality in zip(axes, centralities):
        try:
            pi_HotQCD = np.load("../results/hotqcd/" + centrality + "/pion.npz")
            k_HotQCD = np.load("../results/hotqcd/" + centrality + "/kaon.npz")
            p_HotQCD = np.load("../results/hotqcd/" + centrality + "/proton.npz")
            pi_WB = np.load("../results/wb/" + centrality + "/pion.npz")
            k_WB = np.load("../results/wb/" + centrality + "/kaon.npz")
            p_WB = np.load("../results/wb/" + centrality + "/proton.npz")
        except IOError:
            print("missing file in", system + centrality)
        else:
            pt = pi_HotQCD['pt']
            ratio = pi_WB['v2']/pi_HotQCD['v2']
            err = ratio*np.sqrt((pi_WB['dv2']/pi_WB['v2'])**2 + (pi_HotQCD['dv2']/pi_HotQCD['v2'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Blues(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Blues(0.6), alpha=0.5, linewidth=0)

            pt = k_HotQCD['pt']
            ratio = k_WB['v2']/k_HotQCD['v2']
            err = ratio*np.sqrt((k_WB['dv2']/k_WB['v2'])**2 + (k_HotQCD['dv2']/k_HotQCD['v2'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Oranges(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Oranges(0.6), alpha=0.5, linewidth=0)
           
            pt = p_HotQCD['pt']
            ratio = p_WB['v2']/p_HotQCD['v2']
            err = ratio*np.sqrt((p_WB['dv2']/p_WB['v2'])**2 + (p_HotQCD['dv2']/p_HotQCD['v2'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Greens(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Greens(0.6), alpha=0.5, linewidth=0)

            ax.plot(np.linspace(0,3,100), np.ones(100), linewidth=0.2, color='gray')
            ax.set_ylim([0.85,1.15])

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
    ax1b.set_yticks([0.9, 1.0, 1.1])
    ax3b.annotate('S95/HQCD', xy=(1.025,0.5), va='center', xycoords='axes fraction', rotation=-90, fontsize=texnormal)

    for ax, centrality in zip(axes, centralities):
        try:
            pi_HotQCD = np.load("../results/hotqcd/" + centrality + "/pion.npz")
            k_HotQCD = np.load("../results/hotqcd/" + centrality + "/kaon.npz")
            p_HotQCD = np.load("../results/hotqcd/" + centrality + "/proton.npz")
            pi_s95 = np.load("../results/s95/" + centrality + "/pion.npz")
            k_s95 = np.load("../results/s95/" + centrality + "/kaon.npz")
            p_s95 = np.load("../results/s95/" + centrality + "/proton.npz")
        except IOError:
            print("missing file in", system + centrality)
        else:
            pt = pi_HotQCD['pt']
            ratio = pi_s95['v2']/pi_HotQCD['v2']
            err = ratio*np.sqrt((pi_s95['dv2']/pi_s95['v2'])**2 + (pi_HotQCD['dv2']/pi_HotQCD['v2'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Blues(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Blues(0.6), alpha=0.5, linewidth=0)

            pt = k_HotQCD['pt']
            ratio = k_s95['v2']/k_HotQCD['v2']
            err = ratio*np.sqrt((k_s95['dv2']/k_s95['v2'])**2 + (k_HotQCD['dv2']/k_HotQCD['v2'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Oranges(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Oranges(0.6), alpha=0.5, linewidth=0)
           
            pt = p_HotQCD['pt']
            ratio = p_s95['v2']/p_HotQCD['v2']
            err = ratio*np.sqrt((p_s95['dv2']/p_s95['v2'])**2 + (p_HotQCD['dv2']/p_HotQCD['v2'])**2)
            ax.plot(pt, ratio, linewidth=0.75, color=plt.cm.Greens(0.6))
            ax.fill_between(pt, ratio - err , ratio + err, color=plt.cm.Greens(0.6), alpha=0.5, linewidth=0)
            
            ax.plot(np.linspace(0,3,100), np.ones(100), linewidth=0.2, color='gray')
            
            ax.set_xlabel("$p_T$ [GeV]")
            ax.set_ylim([0.85,1.15])
    
    for ax in fig.axes:
        ax.set_xlim([0,2.5])

    plt.tight_layout(pad=.1, w_pad=0, rect=[0,0,0.975,1])

if __name__== "__main__":
    v2()
