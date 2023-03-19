import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

def angle_sft(x):

    if x < -np.pi/2:
        return x + 2*np.pi
    else:
        return x

def run():

    list_file = [
        'carbonic.csv',
        '../../../H2CO3_CT_H2O_126/330K/carbonic/carbonic.csv',
        '../../../H2CO3_CC_H2O_126/330K/carbonic/carbonic.csv',
    ]

    fig, ax = plt.subplots()

    df_new = None
    for str_file in list_file:
        df_tmp = pd.read_csv(str_file)
        if df_new is None:
            df_new = df_tmp
        else:
            df_new = pd.concat([df_new, df_tmp], ignore_index=True)
    print(df_new)
    df_new = df_new[df_new['dihedral0(rad)'].notnull()]
    alpha = df_new['dihedral0(rad)'].apply(angle_sft)
    beta = df_new['dihedral1(rad)'].apply(angle_sft)

    h, xedges, yedges, image = ax.hist2d(
        alpha,
        beta,
        bins = 100,
        density = True,
        #cmin = 0.001,
        norm = LogNorm(),
        cmap = 'coolwarm'
    )
    fig.colorbar( mappable=image )

    plt.show()

run()
