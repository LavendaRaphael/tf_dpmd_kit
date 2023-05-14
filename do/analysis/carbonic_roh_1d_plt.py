import matplotlib.pyplot as plt
from tf_dpmd_kit import plot
from tf_dpmd_kit import plm
from tf_dpmd_kit import analysis
import pandas as pd
import numpy as np
import matplotlib as mpl

def run(
    ax,
):

    df = pd.read_csv('carbonic_roh_1d.csv', index_col='roh0(ang)')
    for label in ['CC', 'CT', 'TT', 'H2CO3']:
        ser = df[label]
        ser -= min(ser[(ser.index>3.0) & (ser.index<4.2)])
        ax.plot(df.index, df[label], label=label, lw=1)
    
    ax.set_xlabel(r'R($^=$OH$_{proton}$) (Å)')
    ax.set_ylabel('Free energy (kcal/mol)')
    ax.legend()

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots(figsize = (8.6*cm, 5*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_roh_1d',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()
