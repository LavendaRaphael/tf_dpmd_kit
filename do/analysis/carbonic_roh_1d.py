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

    df = analysis.read_multidata([
        '../CC/carbonic/carbonic.product.csv',
        '../CT/carbonic/carbonic.product.csv',
        '../TT/carbonic/carbonic.product.csv',
    ]).dropna()
    print(df)

    bool_c0 = df['dihedral0(rad)'] < np.pi/2
    bool_t0 = df['dihedral0(rad)'] > np.pi/2
    bool_c1 = df['dihedral1(rad)'] < np.pi/2
    bool_t1 = df['dihedral1(rad)'] > np.pi/2
    bool_cc = bool_c0 & bool_c1
    bool_ct = (bool_c0 & bool_t1) | (bool_t0 & bool_c1)
    bool_tt = bool_t0 & bool_t1

    dict_x = {
        'CC': df[bool_cc],
        'CT': df[bool_ct],
        'TT': df[bool_tt],
        'all': df,
    }
    df = pd.DataFrame()
    for label, df_p in dict_x.items():
        np_hist, bin_edges = np.histogram(df_p['roh0(ang)'], bins=200, density=True, range=(0.8,6))
        np_energy = plm.prob_to_deltag(np_hist, temperature=330)
        bin_center = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2
        df_tmp = pd.DataFrame({'roh0(ang)': bin_center, label:np_energy})
        x = df_tmp['roh0(ang)']
        y = df_tmp[label]
        #np_energy -= np.amin(np_energy)
        y -= min(y[(x>3.5) & (x<4.5)])
        ax.plot(x, y, label=label)
        df = pd.concat([df, y], axis=1)
    df = pd.concat([df, x], axis=1)

    ax.set_xlabel(r'R(OH) (Å)')
    ax.set_ylabel('Free energy (kJ/mol)')
    ax.legend()

    #df = pd.DataFrame(data={'roh(ang)': bin_center, 'energy(kJ/mol)': np_energy})
    df.to_csv('carbonic_roh_1d.csv', index=False)

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots(figsize = (8.6*cm, 5*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_roh_1d',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()
