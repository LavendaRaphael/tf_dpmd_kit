import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from tf_dpmd_kit import plm
from tf_dpmd_kit import plot

def angle_sft(x):

    if x < -np.pi/2:
        return x + 2*np.pi
    else:
        return x

def prob2energy(
    prob: float,
    temperature: float, # kelvin
):

    KbT = plm.T2KbT(temperature)
    return -KbT*np.log(prob)

def run(
    list_file: list,
):

    plot.set_rcparam()
    cm = 1/2.54

    fig, ax = plt.subplots()

    df_new = None
    for str_file in list_file:
        print(str_file)
        df_tmp = pd.read_csv(str_file)
        if df_new is None:
            df_new = df_tmp
        else:
            df_new = pd.concat([df_new, df_tmp], ignore_index=True)
    print(df_new)

    df_new = df_new[df_new['dihedral0(rad)'].notnull()]
    alpha = df_new['dihedral0(rad)'].apply(angle_sft)
    beta = df_new['dihedral1(rad)'].apply(angle_sft)

    h, xedges, yedges = np.histogram2d(alpha, beta, bins=100, density=True)
    energy = prob2energy(h, temperature=330)
    energy -= np.amin(energy)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    cmap_new = plt.get_cmap('coolwarm', 10)
    #cmap_new = 'coolwarm'
    image = ax.imshow(
        energy.T, 
        origin = 'lower',
        extent = extent,
        cmap = cmap_new,
        aspect = 'auto',
    )
    fig.colorbar( mappable=image )

    ax.set_xlabel(r'$\alpha$ (rad)')
    ax.set_ylabel(r'$\beta$ (rad)')

    ax.set_xticks([0, np.pi/2, np.pi])
    ax.set_yticks([0, np.pi/2, np.pi])
    ax.set_xticklabels([0, r'$\pi$/2', r'$\pi$'])
    ax.set_yticklabels([0, r'$\pi$/2', r'$\pi$'])

    plot.save(
        fig,
        str_save = 'carbonic_dihedrals_2d',
        tup_size = (8.6*cm, 6.8*cm),
        list_type = ['svg', 'pdf'],
    )

    plt.show()

run(
    list_file = [
        'carbonic.product.csv',
        '../../../H2CO3_CT_H2O_126/330K/carbonic/carbonic.product.csv',
        '../../../H2CO3_CC_H2O_126/330K/carbonic/carbonic.product.csv',
    ]
)
