import matplotlib.pyplot as plt
from tf_dpmd_kit import plot
import pandas as pd
import numpy as np

plot.set_rcparam()
cm = 1/2.54

def run():

    df_data = pd.read_csv('carbonic_intermiccime.csv')
    df_data[df_data > 2] = np.nan

    fig, ax = plt.subplots()

    list_header = [ 'CT','CC', 'TT']
    #list_header = ['CC']
    for header in list_header:
        ax.hist(df_data[header], bins=100, label=header, alpha=0.5, range=(0,2))
        #ax.hist(df_data[header], bins=100, label=header, alpha=0.5, density=True)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Count')
    ax.legend()

    plot.save(
        fig,
        tup_size = (8.6*cm, 5*cm),
        str_save = 'carbonic_intermiccime',
        list_type = ['pdf', 'svg'],
    )

run()

plt.show()

