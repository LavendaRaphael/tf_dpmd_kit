import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tf_dpmd_kit import plm
from tf_dpmd_kit import plot
import matplotlib as mpl

plot.set_rcparam()
cm = 1/2.54
mpl.rcParams['figure.dpi'] = 300

def run(
    fig,
    ax,
):

    df_data = pd.read_csv('carbonic_dihedrals.csv')

    h, xedges, yedges = np.histogram2d(df_data['dihedral0(rad)'], df_data['dihedral1(rad)'], bins=100, density=True)

    energy = plm.prob_to_deltag(h, temperature=330)
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

    #ax.set_aspect(1)
    ax.set_xlabel(r'$\alpha$ (rad)')
    ax.set_ylabel(r'$\beta$ (rad)')

    ax.set_xticks([0, np.pi/2, np.pi])
    ax.set_yticks([0, np.pi/2, np.pi])
    ax.set_xticklabels([0, r'$\pi$/2', r'$\pi$'])
    ax.set_yticklabels([0, r'$\pi$/2', r'$\pi$'])

def main():

    fig, ax = plt.subplots(figsize = (8.6*cm, 7*cm))

    run(fig, ax)

    plot.save(
        fig,
        file_save = 'carbonic_dihedrals_2d',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()



