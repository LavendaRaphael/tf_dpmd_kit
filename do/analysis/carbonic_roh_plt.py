import matplotlib.pyplot as plt
import matplotlib as mpl
from tf_dpmd_kit import analysis
from tf_dpmd_kit import plot
from tf_dpmd_kit import plm
import pandas as pd
import numpy as np
import matplotlib.colors as colors

def run(fig,ax):

    df_data = analysis.read_multidata([
        '../CC/carbonic/carbonic.product.csv',
        '../CT/carbonic/carbonic.product.csv',
        '../TT/carbonic/carbonic.product.csv',
    ]).dropna()

    dh0 = df_data['dihedral0(rad)']
    dh1 = df_data['dihedral1(rad)']

    pio2 = np.pi/2
    #df_data = df_data[ (dh0 > -pio2) & (dh0 < pio2) & (dh1 > -pio2) & (dh1 < pio2) ]
    #df_data = df_data[ ((dh0 > -pio2) & (dh0 < pio2) & ((dh1 < -pio2) | (dh1 > pio2))) | (((dh0 < -pio2) | (dh0 > pio2)) & ((dh1 > -pio2) & (dh1 < pio2))) ]
    print(df_data)
    df_sym = df_data.rename(columns={'roh0(ang)': 'roh1(ang)', 'roh1(ang)': 'roh0(ang)'})
    df_data = pd.concat([df_data, df_sym], ignore_index=True)

    h, xedges, yedges = np.histogram2d(df_data['roh0(ang)'], df_data['roh1(ang)'], bins=200, density=True)
    #h, xedges, yedges = np.histogram2d(df_data['dihedral0(rad)'], df_data['dihedral1(rad)'], bins=200, density=True)

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
    colorbar.ax.set_ylabel('Free energy (kJ/mol)')

    ax.set_xlabel(r'$R$ (Å)')
    ax.set_ylabel(r'$R$ (Å)')

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
