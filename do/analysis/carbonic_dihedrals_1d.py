import matplotlib.pyplot as plt
from tf_dpmd_kit import plot
from tf_dpmd_kit import plm
import pandas as pd
import numpy as np
import matplotlib as mpl

def run(
    ax,
):

    df_data = pd.read_csv('carbonic_dihedrals.csv')
    ser_sum = df_data['dihedral0(rad)'] + df_data['dihedral1(rad)']

    np_hist, bin_edges = np.histogram(ser_sum, bins=100, density=True)
    np_energy = plm.prob_to_deltag(np_hist, temperature=330)
    np_energy -= np.amin(np_energy)
    bin_center = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2

    ax.plot(bin_center, np_energy)

    ax.set_xlabel(r'$\alpha+\beta$ (rad)')
    ax.set_ylabel('Free energy (kJ/mol)')
    ax.set_xticks([0, np.pi/2, np.pi, np.pi*1.5, np.pi*2])
    ax.set_xticklabels([0, r'$\pi$/2', r'$\pi$', r'3$\pi$/2', r'2$\pi$'])
    ax.set_ylim(None, None)

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots(figsize = (8.6*cm, 5*cm))

    run(ax)

    plot.save(
        fig,
        #file_save = 'carbonic_dihedrals_1d',
        list_type = ['pdf', 'svg']
    )

    plt.show()
 

main()