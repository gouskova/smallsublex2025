import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.font_manager
import os
import numpy as np
import pandas as pd
import nclasses as pnc 

def plot_sim_with_ci(values, ci, abline, bins=6, show=True, fname='simulation', ftype='pdf', color=False):
    '''
    plot the maximum size cap in syllables for each simulation
    get confidence intervals for each simulation
    point to where the sublexicon sits within that span
    '''
    libfont = {'fontname':'Linux Libertine O', 'size': 'x-large'}
    sns.set_theme(style='whitegrid')
    h_color='gray'
    fig = sns.histplot(values, discrete=True, color=h_color).set_title("Monte Carlo max size: "+" ".join(fname.split("_")), **libfont)
    if color:
        ci_color='r'
        abline_color='green'
    else:
        ci_color='black'
        abline_color='black'
    plt.axvline(ci[0], color=ci_color)
    plt.axvline(ci[1], color=ci_color)
    plt.axvline(abline, linestyle="--", color=abline_color)
    if show:
        plt.show()
    fig.figure.savefig(os.path.join(os.path.expanduser('~/git/smallsublex/plots'), '.'.join([fname, ftype])))
    plt.clf()
    plt.close()

def plot_syllcounts(fpath, show=True, ftype="pdf", color=True, featpath=""):
    '''
    quick-and-dirty
    '''
    libfont = {'fontname':'Linux Libertine O', 'size': 'x-large'}
    plotdir = os.path.basename(os.path.split(fpath)[0])
    plottitle = ' '.join(plotdir.split("_")).capitalize()
    values = {}
    sns.set_theme(style='whitegrid')
    sns.set_style("ticks")
    vowels = pnc.get_vowels(**{'featpath':featpath})
    wlenths = []
    colors: {}
    with open(fpath, 'r') as f:
        for line in f:
            lenth = len([x for x in line.split() if x in vowels])
            label = f'{str(lenth)}$\sigma$'
            wlenths.append(label)
            if label in values:
                values[label]+=1
            else:
                values[label]=1
    if color:
        colors = {f"{x}$\sigma$": "gray" for x in range(10)}
        colors['0$\sigma$'] = 'blue'
        colors['1$\sigma$'] = 'black'
        colors['2$\sigma$'] = 'yellow'
        colors['3$\sigma$'] = 'green'
        colors['4$\sigma$'] = 'red'
    elif color == 'gr':
        colors='Greys'
    else:
        colors = {f"{x}$\sigma$": "lightgray" for x in range(10)}
        colors['0$\sigma$'] = 'whitesmoke'
        colors['1$\sigma$'] = 'black'
        colors['2$\sigma$'] = 'dimgray'
        colors['3$\sigma$'] = 'darkgray'
        colors['4$\sigma$'] = 'gray'
    xvals = sorted(values)
    yvals = [values[z] for z in xvals]
    fig = sns.catplot(x=xvals, y=yvals, kind='bar', hue=xvals, hue_order=sorted(xvals), legend=False, palette=colors)
    fig.set_axis_labels(plottitle, 'Count', **libfont)
    for c in fig.ax.containers:
            fig.ax.bar_label(c, label_type="edge")
    if show:
        plt.show()
    fig.savefig(os.path.join(os.path.expanduser('~/git/smallsublex/plots'), '.'.join([plotdir, ftype])))
    plt.close()

def plot_syllcount_by_freq(fpath, show=True, ftype="pdf"):
    '''
    input is a tab-separated file, in IPA, with transcriptions space-separated.
    '''
    libfont = {'fontname':'Linux Libertine O', 'size': 'x-large'}
    plotdir = os.path.basename(os.path.split(fpath)[0])
    plottitle = "Length by frequency" 
    values = {}
    sns.set_theme(style='whitegrid')
    colors: {}
    df = pd.read_csv(fpath, sep="\t")
    df['rank']=df.index+1
    print(df.describe())
    fig = sns.relplot(y="syllables", x="rank", data=df)
    if show:
        plt.show()
    fig.figure.savefig(os.path.join(os.path.expanduser('~/git/smallsublex/plots'),"length_by_freq.pdf"))
    plt.close()




if __name__=='__main__':
    #plot_syllcount_by_freq('~/Desktop/nunya_bidnis/LearningData.txt')
    fpath = os.path.expanduser('~/git/smallsublex/data/russian/Features.txt')
    for i in ['freq_noun_stems',
            'freq_adj_stems',
            'astyj_stems_aranea',
            'ost_stems_aranea',
            'freq_astyj',
            'freq_ist',
            'freq_izm',
            'freq_ost',
            'freq_onok',
            'onok_stems_aranea']:
        plot_syllcounts(os.path.join(os.path.expanduser("~/git/smallsublex/data/russian"), i, "LearningData.txt"), color=False, featpath=fpath)
        pass

    fpath = os.path.expanduser('~/git/smallsublex/data/english/Features.txt')
    for i in ['freq_all_adj',
              'freq_en',
              'freq_ify',
              'freq_ize',
              'freq_nouns',
              ]:
        plot_syllcounts(os.path.join(os.path.expanduser("~/git/smallsublex/data/english"), i, "LearningData.txt"),
                        color=False, featpath=fpath)
