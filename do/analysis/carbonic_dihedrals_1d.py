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

def run():

    list_file = [
        'carbonic.csv',
        '../../../H2CO3_CT_H2O_126/330K/carbonic/carbonic.csv',
        '../../../H2CO3_CC_H2O_126/330K/carbonic/carbonic.csv',
    ]

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
    ser_alpha = df_new['dihedral0(rad)'].apply(angle_sft)
    ser_beta = df_new['dihedral1(rad)'].apply(angle_sft)
    ser_sum = ser_alpha + ser_beta

    np_hist, bin_edges = np.histogram(ser_sum, bins=100, density=True)
    np_energy = prob2energy(np_hist, temperature=330)
    np_energy -= np.amin(np_energy)
    bin_center = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2

    plt.plot(bin_center, np_energy, color='tab:grey')
 
    ax.set_xlabel(r'$\alpha+\beta$ (rad)')
    ax.set_ylabel('Free energy (kJ/mol)')
    ax.set_xticks([0, np.pi/2, np.pi, np.pi*1.5, np.pi*2])
    ax.set_xticklabels([0, r'$\pi$/2', r'$\pi$', r'3$\pi$/2', r'2$\pi$'])

    ax.set_ylim(None, 30)

    plot.save(
        fig,
        str_save = 'carbonic_dihedrals_1d',
        tup_size = (8.6*cm, 5*cm),
        list_type = ['svg', 'pdf'],
    )

    plt.show()

    df_save = pd.DataFrame( {'alpha+beta(rad)': bin_center, 'energy(kJ/mol)': np_energy} )
    df_save.to_csv('carbonic_dihedrals_1d.csv', index=False)

run()
