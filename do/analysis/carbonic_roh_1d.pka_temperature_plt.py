import matplotlib.pyplot as plt
from tf_dpmd_kit import plot
from tf_dpmd_kit import plm
from tf_dpmd_kit import analysis
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib.lines import Line2D

def run(
    ax,
):

    list_marker = Line2D.filled_markers[1:]
    mx = 0
    df = pd.read_csv('carbonic_roh_1d.pka.csv')
    for state in ['CC', 'CT', 'TT', 'H2CO3']:
        df_new = df.groupby('state').get_group(state)
        ax.errorbar(df_new['temperature(K)'], df_new['pka'], ls=':', marker=list_marker[mx], capsize=2, lw=1, markersize=2, label=state)
    
    ax.set_xlabel(r'Temperature (K)')
    ax.set_ylabel('pKa')
    ax.legend()

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots(figsize = (8.6*cm, 5*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_roh_1d.pka_temperature',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()
