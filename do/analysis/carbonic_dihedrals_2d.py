import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tf_dpmd_kit import plm
from tf_dpmd_kit import plot
from tf_dpmd_kit import analysis
import matplotlib as mpl
import matplotlib.colors as colors

def run(
    fig,
    ax,
):

    df = analysis.read_multidata([
        '../CC/carbonic/carbonic.product.csv',
        '../CT/carbonic/carbonic.product.csv',
        '../TT/carbonic/carbonic.product.csv',
    ]).dropna()
    df = df[df['ncarbonyl'] == 1.0]

    df_sym = df.rename(columns={'dihedral1(rad)': 'dihedral0(rad)', 'dihedral0(rad)': 'dihedral1(rad)'})
    df = pd.concat([df, df_sym], ignore_index=True)

    bool_c0 = df['dihedral0(rad)'] < np.pi/2
    bool_t0 = df['dihedral0(rad)'] > np.pi/2
    bool_c1 = df['dihedral1(rad)'] < np.pi/2
    bool_t1 = df['dihedral1(rad)'] > np.pi/2
    bool_cc = bool_c0 & bool_c1
    bool_ct = (bool_c0 & bool_t1) | (bool_t0 & bool_c1)
    bool_tt = bool_t0 & bool_t1

    df_p = df[bool_tt]
    h, xedges, yedges = np.histogram2d(df_p['dihedral0(rad)'], df_p['dihedral1(rad)'], bins=200, density=True)

    energy = plm.prob_to_deltag(h, temperature=330)
    energy -= np.amin(energy)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    image = ax.imshow(
        energy.T,
        origin = 'lower',
        extent = extent,
        #cmap = 'Blues_r',
        cmap = 'coolwarm',
        aspect = 'auto',
        #norm =colors.BoundaryNorm(boundaries=np.linspace(0, 30, 11), ncolors=256, extend='max') ,
    )
    colorbar = fig.colorbar( mappable=image, ax=ax, extend='max')
    colorbar.ax.set_ylabel('Free energy (kcal/mol)')

    ax.set_xlabel(r'$\alpha$ (rad)')
    ax.set_ylabel(r'$\beta$ (rad)')

    ax.set_xticks([0, np.pi/2, np.pi])
    ax.set_yticks([0, np.pi/2, np.pi])
    ax.set_xticklabels([0, r'$\pi$/2', r'$\pi$'])
    ax.set_yticklabels([0, r'$\pi$/2', r'$\pi$'])

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots(figsize = (8.6*cm, 7*cm))

    run(fig, ax)

    plot.save(
        fig,
        #file_save = 'carbonic_dihedrals_2d',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()



