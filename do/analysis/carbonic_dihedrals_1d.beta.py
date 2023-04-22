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

    df_data = analysis.read_multidata([
        '../CC/carbonic/carbonic_dihedrals.csv',
        '../CT/carbonic/carbonic_dihedrals.csv',
        '../TT/carbonic/carbonic_dihedrals.csv',
    ])

    df_sym = df_data.rename(columns={'dihedral1(rad)': 'dihedral0(rad)', 'dihedral0(rad)': 'dihedral1(rad)'})
    df_data = pd.concat([df_data, df_sym], ignore_index=True)

    df_data = df_data[df_data['dihedral0(rad)'] > np.pi/2]
    print(df_data)
    ser_sum = df_data['dihedral1(rad)']

    np_hist, bin_edges = np.histogram(ser_sum, bins=100, density=True)
    np_energy = plm.prob_to_deltag(np_hist, temperature=330)
    np_energy -= np.amin(np_energy)
    bin_center = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2

    ax.plot(bin_center, np_energy)

    ax.set_xlabel(r'$\beta [\alpha>\pi/2]$ (rad)')
    ax.set_ylabel('Free energy (kJ/mol)')
    ax.set_xticks([ 0, np.pi/2, np.pi])
    ax.set_xticklabels([ 0, r'$\pi$/2', r'$\pi$'])
    ax.set_xlim(-np.pi/2, np.pi*1.5)
    ax.set_ylim(None, None)

    df = pd.DataFrame(data={'dihedral(rad/pi)': bin_center/np.pi, 'energy(kJ/mol)': np_energy})
    df.to_csv('carbonic_dihedrals_1d.beta.csv', index=False)

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots(figsize = (8.6*cm, 3*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_dihedrals_1d.beta',
        list_type = ['pdf', 'svg']
    )

    plt.show()
 

main()
