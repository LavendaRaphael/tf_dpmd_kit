import matplotlib.pyplot as plt
import matplotlib as mpl
from tf_dpmd_kit import plot
from tf_dpmd_kit import plm
import os
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D


def run(
    ax1, ax2,
):

    dict_color = {
        'CC': 'tab:blue',
        'CT': 'tab:orange',
        'TT': 'tab:green',
        'HCO3': 'tab:purple',
    }
    dict_marker = {
        'CC': 'o',
        'CT': 'v',
        'TT': '^',
        'HCO3': '>',
    }

    list_header = ['CC', 'CT', 'TT','HCO3']

    file_data = 'carbonic_statistic.temperature.csv'
    dfgb = pd.read_csv(file_data, index_col=['state']).groupby(level='state')
    ser_temperature = dfgb.get_group(list_header[0])['temperature(K)']
    for header in list_header:
        color = dict_color[header]
        marker = dict_marker[header]
        df = dfgb.get_group(header)
        ax1.errorbar(ser_temperature, df['prop'], yerr=df['prop_sem'], ls=':', marker=marker, markersize=2, lw=1, color=color, capsize=2)
        ax2.errorbar(ser_temperature, df['prop'], yerr=df['prop_sem'], ls=':', marker=marker, markersize=2, lw=1, color=color, capsize=2)

    ax1.set_ylim(0.75,None)
    ax2.set_ylim(None,0.25)
    ax1.spines.bottom.set_visible(False)
    ax1.xaxis.set_ticks_position('none')
    ax2.spines.top.set_visible(False)

    kwargs = dict(
        marker=[(-1, -.5), (1, .5)], 
        markersize=6,
        linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

    ax2.set_xlim(None, 375)
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Proportion', y= 1.)

    plot.add_text(
        ax1,
        dict_text = {
            (355, 0.8): 'CC',
        }
    )
    plot.add_text(
        ax2,
        dict_text = {
            (355,  0.1): 'CT',
            (355,  0.05): r'HCO$_3^-$',
            (355,  0.0): 'TT',
        }
    )

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300
    mpl.rcParams['figure.constrained_layout.use'] = False

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize = (4.3*cm, 4*cm), gridspec_kw={'left':0.23, 'right': 0.95, 'bottom':0.2, 'top':0.95})
    #fig.subplots_adjust(hspace=0.005)

    run(ax1, ax2)

    plot.save(
        fig,
        file_save = 'carbonic_statistic.temperature_prop',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()


