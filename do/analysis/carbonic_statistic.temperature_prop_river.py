import matplotlib.pyplot as plt
import matplotlib as mpl
from tf_dpmd_kit import plot
from tf_dpmd_kit import plm
import os
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D

plot.set_rcparam()
cm = 1/2.54
mpl.rcParams['figure.dpi'] = 300

def run(
    ax,
):

    dict_color = {
        'TT': 'tab:blue',
        'CT': 'tab:orange',
        'CC': 'tab:green',
        'H2CO3': 'tab:red',
        'HCO3': 'tab:purple',
    }

    file_data = 'carbonic_statistic.temperature.csv'
    list_header = ['TT', 'CT', 'CC', 'HCO3']

    dfgb = pd.read_csv(file_data, index_col=['state']).groupby(level='state')
    ser_temperature = dfgb.get_group(list_header[0])['temperature(K)']
    ser_0 = pd.Series([0]*len(ser_temperature))
    for header in list_header:
        ser_1 = dfgb.get_group(header)['prop'].reset_index(drop=True)
        if header in dict_color:
            color = dict_color[header]
        ser_1 = ser_0 + ser_1
        ax.fill_between(ser_temperature, ser_0, ser_1, lw=1, color=color, alpha=0.5)
        ser_0 = ser_1

    plot.add_text(
        ax,
        dict_text = {
            (330, 0.65): 'TT',
            (330, 0.87): 'CT',
            (340, 0.87): 'CC',
            (330, 0.97): r'HCO$_3^-$',
        },
        va = 'center',
        ha = 'center',
        color = 'white',
        fontweight = 'bold',
        fontsize = 8,
    ),
    plot.add_arrow(
        ax,
        dict_arrow = {
            'x': [(342, 0.87), (348, 0.92)]
        },
        arrowstyle = 'simple, head_length=2, head_width=2, tail_width=0.1',
        color = 'white'
    )

    ax.set_xlim(290, 350)
    ax.set_ylim(0.5, 1)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Proportion'),

def main():

    fig, ax = plt.subplots( figsize = (8.6*cm, 5*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_statistic.temperature_prop_river',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()


