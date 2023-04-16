import matplotlib.pyplot as plt
import matplotlib as mpl
from tf_dpmd_kit import plot
from tf_dpmd_kit import plm
import os
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D


def run(
    ax,
):

    dict_color = {
        'TT': 'tab:blue',
        'CT': 'tab:orange',
        'CC': 'tab:green',
        'HCO3': 'tab:purple',
    }
    dict_marker = {
        'TT': 'o',
        'CT': 'v',
        'CC': '^',
        'HCO3': '>',
    }

    list_header = ['TT', 'CT', 'CC', 'HCO3']

    file_data = 'carbonic_statistic.temperature.csv'
    dfgb = pd.read_csv(file_data, index_col=['state']).groupby(level='state')
    ser_temperature = dfgb.get_group(list_header[0])['temperature(K)']
    for header in list_header:
        color = dict_color[header]
        marker = dict_marker[header]
        df = dfgb.get_group(header)
        ax.errorbar(ser_temperature, df['freqprop'], yerr = df['freqprop_sem'], ls=':', marker=marker, markersize=2, lw=1, color=color, capsize=2)

    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Frequency proportion')

    plot.add_text(
        ax,
        dict_text = {
            (355, 0.2): 'TT',
            (355, 0.3): 'CT',
            (355, 0.4): r'HCO$_3^-$',
            (355, 0.08): 'CC',
        }
    )
    ax.set_xlim(None, 375)

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots( figsize = (4.3*cm, 4*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_statistic.temperature_freqprop',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()


