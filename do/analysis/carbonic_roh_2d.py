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

    print(df_data)
    #df_sym = df_data.rename(columns={'roh0(ang)': 'roh1(ang)', 'roh1(ang)': 'roh0(ang)'})
    #df_data = pd.concat([df_data, df_sym], ignore_index=True)

    h, xedges, yedges = np.histogram2d(df_data['roh0(ang)'], df_data['roh1(ang)'], bins=[300, 300], density=True, range=[[0.8, 10],[0.8, 1.3]])

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

    ax.set_xlabel(r'R($^=$OH$_{proton}$) (Å)')
    ax.set_ylabel(r'R(O$_{OH}$H$_{OH}$) (Å)')

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots(figsize = (8.6*cm, 5*cm))

    run(fig, ax)

    plot.save(
        fig,
        file_save = 'carbonic_roh_2d',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()
