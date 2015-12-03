#!/usr/bin/env python3

import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import glob
import random


aspect = 1/1.618
resolution = 72.27
columnwidth = 246/resolution
textwidth = 510/resolution
textiny, texsmall, texnormal = 8.0, 9.25, 10.0
offblack = '#262626'

plt.rcdefaults()
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['CMU Serif'],
    'font.size': texsmall,
    'legend.fontsize': texsmall,
    'axes.labelsize': texsmall,
    'axes.titlesize': texsmall,
    'xtick.labelsize': textiny,
    'ytick.labelsize': textiny,
    'font.weight': 400,
    'axes.labelweight': 400,
    'axes.titleweight': 400,
    'lines.linewidth': .9,
    'lines.markersize': 3,
    'lines.markeredgewidth': .1,
    'patch.linewidth': .9,
    'axes.linewidth': .5,
    'xtick.major.width': .5,
    'ytick.major.width': .5,
    'xtick.minor.width': .5,
    'ytick.minor.width': .5,
    'xtick.major.size': 2,
    'ytick.major.size': 2,
    'xtick.minor.size': 1.3,
    'ytick.minor.size': 1.3,
    'xtick.major.pad': 1.8,
    'ytick.major.pad': 1.8,
    'text.color': offblack,
    'axes.edgecolor': offblack,
    'axes.labelcolor': offblack,
    'xtick.color': offblack,
    'ytick.color': offblack,
    'legend.numpoints': 1,
    'legend.scatterpoints': 1,
    'legend.frameon': False,
    'image.interpolation': 'none',
    'pdf.fonttype': 42
})


plot_functions = {}


def plot(f):
    def wrapper(*args, **kwargs):
        print(f.__name__)
        f(*args, **kwargs)
        plt.savefig('{}.pdf'.format(f.__name__))
        plt.close()

    plot_functions[f.__name__] = wrapper

    return wrapper


def finish(despine=True, remove_ticks=False, pad=0.1, h_pad=None, w_pad=None):
    fig = plt.gcf()

    for ax in fig.axes:
        if despine:
            for spine in 'top', 'right':
                ax.spines[spine].set_visible(False)
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')
        else:
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

    fig.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)


def set_loc(ax, xy=None, nbins=5, steps=[1, 2, 3, 4, 10],
            prune=None, minor=0):
    if xy == 'x':
        axes = ax.xaxis,
    elif xy == 'y':
        axes = ax.yaxis,
    else:
        axes = ax.xaxis, ax.yaxis

    for axis in axes:
        axis.set_major_locator(
            tick.MaxNLocator(nbins=nbins, steps=steps, prune=prune)
        )
        if minor:
            axis.set_minor_locator(tick.AutoMinorLocator(minor))


@plot
def eos_compare():

    # unit conversion
    def max_trace(thermal):
        e, p, s, T = thermal.T
        return np.max((e - 3*p) / T**4 * 0.19733**3)
    
    # figure properties
    fig, axes = plt.subplots(nrows=4, figsize=(columnwidth, 3*aspect*columnwidth), sharex=True)
    
    # EoS labels
    colors = plt.cm.Blues(0.6), plt.cm.Reds(0.6), plt.cm.Greens(0.6)
    shapes = 's', '^', '*'
    sizes = 4, 5, 7
    labels = 'HQCD', 'WB', 'S95'
    shifts = 0, 5, -30

    # plot mean pt
    charged = np.sort([np.load(eos)['charged'] for eos in glob.glob('mean-pt/results/samples/*')])
    charged_err = np.sort([np.load(eos)['charged_err'] for eos in glob.glob('mean-pt/results/samples/*')])
    axes[0].errorbar(np.arange(len(charged)), charged, yerr=2*charged_err, fmt='o', color=plt.cm.Blues(0.6), capsize=0, lw=1, ms=4, mew=0.2, mec='w')
    perc = np.percentile(charged, [2.5,97.5])
    axes[0].fill_between(np.linspace(0,100,10), perc[0]*np.ones(10), perc[1]*np.ones(10), color='0.05', alpha=0.2, lw=0)
    axes[0].annotate(r'$\pm 2 \sigma$', xy=(97.5, 0.5755), va='top', ha='right', color='0.4')
    axes[0].annotate("$\pm{:.1f}$%".format(round(100*(perc[1]-perc[0])/(perc[1]+perc[0]),2)), xy=(50, perc[1]), ha='center', va='bottom')

    for eos, offset, color, shape, size, label, justify, shift in zip(['hotqcd', 'wb', 's95'], [120,140,160], colors, shapes, sizes, labels, ['left', 'left', 'right'], [4,4,-4]):
        charged = np.load('mean-pt/results/' + eos + '.npz')['charged']
        charged_err = np.load('mean-pt/results/' + eos + '.npz')['charged_err']
        axes[0].errorbar(offset, charged, yerr=2*charged_err, fmt=shape, color=color, capsize=0, lw=1, ms=size, mew=0)
        axes[0].set_ylabel(r'$p_T$ [GeV]')
        HQCD_pt = charged if eos == 'hotqcd' else HQCD_pt
        ratio = 100*(charged/HQCD_pt - 1) if eos != 'hotqcd' else 0
        if eos != 'hotqcd':
            axes[0].annotate("${:+.1f}$%".format(round(ratio,2)), xy=(offset+shift, charged), ha=justify, va='center')
    
    axes[0].set_ylim(0.595, 0.625)
    axes[0].set_yticks([0.6, 0.61, 0.62])
    
    # plot v2, v3
    HQCD_v2 = 0
    HQCD_v3 = 0
    for ax, flowfile, flowlabel in zip(axes[1:3], ['.v2.dat', '.v3.dat'], [r'$v_2\{2\}$', r'$v_3\{2\}$']):
        
        charged = np.array([np.load('mean-pt/results/samples/hotqcd_{}'.format(ieos) + '.npz')['charged'] for ieos in range(0,100)])
        vn  = np.array([np.loadtxt('flow/results/hotqcd_{}'.format(ieos) + flowfile) for ieos in range(0,100)])
        vn = vn[np.argsort(np.column_stack((charged, vn))[:,0])]
        ax.errorbar(np.arange(len(vn)), vn[:,0], yerr=2*vn[:,1], fmt='o', color=plt.cm.Blues(0.6), capsize=0, lw=1, ms=4, mew=0.2, mec='w', label='HQCD samples')
        perc = np.percentile(vn[:,0], [2.5,97.5])
        ax.fill_between(np.linspace(0,100,10), perc[0]*np.ones(10), perc[1]*np.ones(10), color='0.05', alpha=0.2, lw=0)
        ax.annotate("$\pm{:.1f}$%".format(round(100*(perc[1]-perc[0])/(perc[1]+perc[0]),2)), xy=(50, perc[1]), ha='center', va='bottom')

        for eos, offset, color, shape, size, label, justify, shift, label in \
            zip(['hotqcd', 'wb', 's95'], [120, 140, 160], colors, shapes, sizes, labels, ['left', 'left', 'right'], [4,4,-4], ['HQCD', 'WB', 'S95']):
            vn = np.array(np.loadtxt('flow/results/' + eos + flowfile))
            ax.errorbar(offset, vn[0], yerr=2*vn[1], fmt=shape, color=color, capsize=0, lw=1, ms=size, mew=0, label=label)
            ax.set_ylabel(flowlabel)
            HQCD_vn = vn[0] if eos == 'hotqcd' else HQCD_vn
            ratio = 100*(vn[0]/HQCD_vn - 1) if eos != 'hotqcd' else 0
            if eos != 'hotqcd':
                ax.annotate("${:+.1f}$%".format(round(ratio,2)), xy=(offset+shift, vn[0]), ha=justify, va='center')

    # plot max I(T)
    charged = np.array([np.load('mean-pt/results/samples/hotqcd_{}'.format(ieos) + '.npz')['charged'] for ieos in range(100)])
    max = np.array([max_trace(np.loadtxt('../splines/eos-tables/hotqcd_{}.dat'.format(ieos))) for ieos in range(100)])
    max = max[np.argsort(np.column_stack((charged, max))[:,0])]
    perc = np.percentile(max, [2.5,97.5])
    axes[3].errorbar(np.arange(100), max, fmt='o', color=plt.cm.Blues(0.6), capsize=0, lw=1, ms=4, mew=0.2, mec='w')
    axes[3].fill_between(np.linspace(0,100,10), perc[0]*np.ones(10), perc[1]*np.ones(10), color='0.05', alpha=0.2, lw=0)
    axes[3].annotate("$\pm{:.1f}$%".format(round(100*(perc[1]-perc[0])/(perc[1]+perc[0]),2)), xy=(50, perc[1]), ha='center', va='bottom')
    
    for eos, offset, color, shape, size, label, justify, shift in zip(['hotqcd', 'wb', 's95'], [120, 140, 160], colors, shapes, sizes, labels, ['left', 'left', 'right'], [4,4,-4]):
        max  = max_trace(np.loadtxt('../eos-from-measure/' + eos + '.dat'))
        axes[3].errorbar(offset, max, fmt=shape, color=color, capsize=0, lw=1, ms=size, mew=0, label=label)
        HQCD_I = max if eos == 'hotqcd' else HQCD_I
        ratio = 100*(max/HQCD_I- 1) if eos != 'hotqcd' else 0
        if eos != 'hotqcd':
            axes[3].annotate("${:+.1f}$%".format(round(ratio,2)), xy=(offset+shift, max), ha=justify, va='center')
            
    axes[3].annotate('200 GeV Au+Au\n20–30% centrality', xy=(0,5.5), ha='left', va='center')

    axes[1].set_ylim(0.051,0.059)
    axes[1].set_yticks([0.052, 0.054, 0.056, 0.058])   
    axes[2].set_ylim(0.017,0.023)
    axes[2].set_yticks([0.018, 0.020, 0.022])
    axes[3].set_xlim(-10,170)
    axes[3].set_ylim(3,6)
    axes[3].set_xticks(np.arange(0,180,20))
    axes[3].set_xticklabels(['0', '20', '40', '60', '80', '100', 'HQCD', 'WB', 'S95'], rotation=45)
    axes[3].set_xlabel('EoS Number', labelpad=0)
    axes[3].set_ylabel(r'max $I(T)$')
    axes[3].set_yticks([3.5, 4.5, 5.5])
    
    handles, labels = axes[2].get_legend_handles_labels()
    axes[0].legend(handles[::-1], labels[::-1], loc=3)

    finish()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('plots', nargs='*')
    args = parser.parse_args()

    if args.plots:
        for i in args.plots:
            if i.endswith('.pdf'):
                i = i[:-4]
            if i in plot_functions:
                plot_functions[i]()
            else:
                print('unknown plot:', i)
    else:
        for f in plot_functions.values():
            f()


if __name__ == "__main__":
    main()
