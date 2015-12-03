#!/usr/bin/env python


from __future__ import division, print_function, unicode_literals

import functools

import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib import cm
import numpy as np
import os, sys
import h5py
import scipy
from scipy.optimize import leastsq

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
def hbt():
    
    blue = '#3399FF'
    magenta = '#B8008A'
    orange = '#FF9900'

    row_labels = '$R_s$ [fm]', '$R_o$ [fm]', '$R_l$ [fm]'
    HBT_radii = 'Ro', 'Rs', 'Rl'
    HBT_errors = 'dRo', 'dRs', 'dRl'
    HBT_indices = 4, 7, 10
    column_labels = '0-10%', '10-20%', '20-40%'
    centralities = '/0to10/', '/10to20/', '/20to40/'
    experimental = '/0to10_pi+.dat', '/10to20_pi+.dat', '/20to40_pi+.dat'
    mass = 0.1396

    fig = plt.figure(figsize=(textwidth, textwidth*aspect))
    gs = plt.GridSpec(4, 9, wspace=0.0, hspace=0.0)

    for row, row_label, R, dR, HBT_index in zip(range(4), row_labels, HBT_radii, HBT_errors, HBT_indices):
        for col, col_label, centrality, exp in zip(range(3), column_labels, centralities, experimental):
            
            # initialize subfigure environment
            ax = fig.add_subplot(gs[row,col*3:(col+1)*3])

            # load model data
            try:
                R_HQCD = np.load("../results-staged/hotqcd/" + centrality + "/R.npz")
                R_WB = np.load("../results-staged/wb/" + centrality + "/R.npz")
                R_S95 = np.load("../results-staged/s95/" + centrality + "/R.npz")
                R_exp = np.loadtxt("phenix-hbt" + exp)
            except IOError:
                print("missing file")
            else:
                mt = np.sqrt(R_HQCD['kt_list']**2 + mass**2)
                R_ = R_HQCD[R]
                print(R_)
                dR_ = R_HQCD[dR]
                ax.plot(mt, R_, linewidth=0.75, color=magenta, label='HQCD $\pi^\pm$', zorder=2)
                ax.fill_between(mt, R_ - 2*dR_, R_ + 2*dR_, alpha=0.3, linewidth=0, color=magenta, zorder=2)
                
                mt = np.sqrt(R_WB['kt_list']**2 + mass**2)
                R_ = R_WB[R]
                dR_ = R_WB[dR]
                ax.plot(mt, R_, linewidth=0.75, color=orange, label='WB $\pi^\pm$', zorder=3)
                ax.fill_between(mt, R_ - 2*dR_, R_ + 2*dR_, alpha=0.3, linewidth=0, color=orange, zorder=3)
                
                mt = np.sqrt(R_S95['kt_list']**2 + mass**2)
                R_ = R_S95[R]
                dR_ = R_S95[dR]
                ax.plot(mt, R_, linewidth=0.75, color=blue, label='S95 $\pi^\pm$', zorder=1)
                ax.fill_between(mt, R_ - 2*dR_, R_ + 2*dR_, alpha=0.3, linewidth=0, color=blue, zorder=1)
                
                mt = R_exp[:,0]
                R_ = R_exp[:,HBT_index]
                dR_ = R_*(R_exp[:,HBT_index + 2]/100.)
                ax.errorbar(mt, R_, yerr=2.*dR_, fmt='o', markersize=3, color=offblack, markerfacecolor=offblack, markeredgecolor='none', capsize=0, zorder=4, label='PHENIX $\pi^+$')
          
            # x-axis
            ax.set_xlim([0.0, 1.25])

            # y-axis
            ax.set_ylim([0, 10])
            if(col == 0):
                ax.set_yticks(np.arange(0, 10, 2))
                ax.set_ylabel(row_label)
            if(row == 0):
                ax.annotate(col_label, xy=(0.5,0.9), xycoords='axes fraction', ha='center', fontsize=texnormal)
            
            # all axes
            despine(ax)
            if(row != 3):
                ax.tick_params(labelbottom='off')
            if(col != 0):
                ax.tick_params(labelleft='off')
            if(row == 0 and col==2): 
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.18,1.05), handletextpad=0.5, handlelength=1.4, labelspacing=0)
        
    # calculate R_out / R_side
    experimental = '/0to10_pi+_ratio.dat', '/10to20_pi+_ratio.dat', '/20to40_pi+_ratio.dat'
    for col, centrality, exp in zip(range(3), centralities, experimental):    
        
        # initialize subfigure environment
        ax = fig.add_subplot(gs[3,col*3:(col+1)*3])

        # load model data
        try:
            R_HQCD = np.load("../results-staged/hotqcd/" + centrality + "/R.npz")
            R_WB = np.load("../results-staged/wb/" + centrality + "/R.npz")
            R_S95 = np.load("../results-staged/s95/" + centrality + "/R.npz")
            R_exp = np.loadtxt("phenix-hbt" + exp)
        except IOError:
            print("missing file")
        else:
            mt = np.sqrt(R_HQCD['kt_list']**2 + mass**2)
            ratio = R_HQCD['Rs']/R_HQCD['Ro']
            error = ratio*np.sqrt((R_HQCD['dRs']/R_HQCD['Rs'])**2 + (R_HQCD['dRo']/R_HQCD['Ro'])**2)
            ax.plot(mt, ratio, linewidth=0.75, color=magenta, label='HQCD', zorder=2)
            ax.fill_between(mt, ratio - 2*error, ratio + 2*error, linewidth=0, color=magenta, alpha=0.3, zorder=2)

            mt = np.sqrt(R_WB['kt_list']**2 + mass**2)
            ratio = R_WB['Rs']/R_WB['Ro']
            error = ratio*np.sqrt((R_WB['dRs']/R_WB['Rs'])**2 + (R_WB['dRo']/R_WB['Ro'])**2)
            ax.plot(mt, ratio, linewidth=0.75, color=orange, label='WB', zorder=3)
            ax.fill_between(mt, ratio - 2*error, ratio + 2*error, linewidth=0, color=orange, alpha=0.3, zorder=3)
            
            mt = np.sqrt(R_S95['kt_list']**2 + mass**2)
            ratio = R_S95['Rs']/R_S95['Ro']
            error = ratio*np.sqrt((R_S95['dRs']/R_S95['Rs'])**2 + (R_S95['dRo']/R_S95['Ro'])**2)
            ax.plot(mt, ratio, linewidth=0.75, color=blue, label='S95', zorder=1)
            ax.fill_between(mt, ratio - 2*error, ratio + 2*error, linewidth=0, color=blue, alpha=0.3, zorder=1)
            
            ax.plot(np.linspace(0,1.25,100), np.ones(100), linewidth=0.2, color='gray')

            mt = R_exp[:,0]
            ratio = R_exp[:,1]
            dratio = R_exp[:,2]
            ax.errorbar(mt, ratio, yerr=2*dratio, label='PHENIX', fmt='o', markersize=3, color=offblack, markerfacecolor=offblack, markeredgecolor='none', capsize=0, zorder=3)

        # x-axis
        ax.set_xlim([0.0, 1.25])
        ax.set_xticks(np.arange(0.0, 1.25, 0.25))
        if(col == 2):
            ax.set_xticks(np.arange(0.0, 1.5, 0.25))
        ax.set_xlabel(r"$m_T$ [GeV]")

        # y-axis
        ax.set_ylim([0.6, 1.4])
        if(col == 0):
            ax.set_yticks([0.8, 1.0, 1.2])
            ax.set_ylabel(r"$R_o / R_s$")
        
        # all axes
        despine(ax)
        if(col != 0):
            ax.tick_params(labelleft='off')

    plt.tight_layout(pad=.1, w_pad=0, rect=[0,0,0.975,1])

if __name__== "__main__":
    hbt()
