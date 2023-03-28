import matplotlib.pyplot as plt
from tf_dpmd_kit import plot
from tf_dpmd_kit import plm
import pandas as pd
import numpy as np

def main():

    plot.set_rcparam()
    cm = 1/2.54

    fig, ax = plt.subplots(figsize=(8.6*cm, 5*cm))

    dict_file = {
        '290K': '../290K/carbonic/carbonic_dihedrals.csv',
        '310K': '../310K/carbonic/carbonic_dihedrals.csv',
        '330K': '../330K/carbonic/carbonic_dihedrals.csv',
        '350K': '../350K/carbonic/carbonic_dihedrals.csv',
    }
    dict_temperature = {
        '290K': 290,
        '310K': 310,
        '330K': 330,
        '350K': 350,
    }
    for key, file_data in dict_file.items():
        temperature = dict_temperature[key]
        print(key, file_data, temperature)

        df_data = pd.read_csv(file_data)
        ser_sum = df_data['dihedral0(rad)'] + df_data['dihedral1(rad)']

        np_hist, bin_edges = np.histogram(ser_sum, bins=100, density=True)
        np_energy = plm.prob_to_deltag(np_hist, temperature=temperature)
        np_energy -= np.amin(np_energy)
        bin_center = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2

        ax.plot(bin_center, np_energy, label=key)

    ax.set_xlabel(r'$\alpha+\beta$ (rad)')
    ax.set_ylabel('Free energy (kJ/mol)')
    ax.set_xticks([0, np.pi/2, np.pi, np.pi*1.5, np.pi*2])
    ax.set_xticklabels([0, r'$\pi$/2', r'$\pi$', r'3$\pi$/2', r'2$\pi$'])

    ax.legend()

    plt.show()

main()
