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

def run_(
    ax,
):

    file = 'carbonic_dihedrals.csv'

    temperature=330

    df_energy = pd.DataFrame()
    df_data = pd.read_csv(file)
    ser_sum = df_data['dihedral0(rad)'] + df_data['dihedral1(rad)']
    np_hist, bin_edges = np.histogram(ser_sum, bins=100, range=(-np.pi, 3*np.pi))

    np_energy = plm.prob_to_deltag(np_hist, temperature=temperature)
    np_lower = plm.prob_to_deltag(np_hist-np.sqrt(np_hist), temperature=temperature)
    np_upper = plm.prob_to_deltag(np_hist+np.sqrt(np_hist), temperature=temperature)

    e_min = np.amin(np_energy)
    np_energy -= e_min
    np_lower -= e_min
    np_upper -= e_min

    bin_center = (bin_edges[1:] +  bin_edges[:-1])/2
    ax.plot(bin_center, np_energy)
    ax.fill_between(bin_center, np_lower, np_upper, alpha=0.5)

    ax.set_xlabel(r'$\alpha+\beta$ (rad)')
    ax.set_ylabel('Free energy (kJ/mol)')
    ax.set_xticks([0, np.pi/2, np.pi, np.pi*1.5, np.pi*2])
    ax.set_xticklabels([0, r'$\pi$/2', r'$\pi$', r'3$\pi$/2', r'2$\pi$'])
    ax.set_ylim(None, None)


def run__(
    ax,
):

    list_file = [
        '../TT/carbonic/carbonic_dihedrals.csv',
        '../CT/carbonic/carbonic_dihedrals.csv',
        '../CC/carbonic/carbonic_dihedrals.csv',
    ]

    temperature=330

    df_energy = pd.DataFrame()
    for file in list_file:
        df_data = pd.read_csv(file)
        ser_sum = df_data['dihedral0(rad)'] + df_data['dihedral1(rad)']
        np_hist, bin_edges = np.histogram(ser_sum, bins=100, density=True, range=(-np.pi, 3*np.pi))
        np_energy = plm.prob_to_deltag(np_hist, temperature=temperature)
        np_energy -= np.amin(np_energy)
        ser_energy = pd.Series(np_energy)
        df_energy = pd.concat([df_energy, ser_energy], axis=1)

    df_mean = df_energy.mean(axis=1)
    df_sem = df_energy.std(axis=1)

    bin_center = (bin_edges[1:] +  bin_edges[:-1])/2

    ax.plot(bin_center, df_mean)
    ax.fill_between(bin_center, df_mean-df_sem, df_mean+df_sem, alpha=0.5)

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
