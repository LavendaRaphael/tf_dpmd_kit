from tf_dpmd_kit import plm
from tf_dpmd_kit import analysis
import pandas as pd
import numpy as np

def run(
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
        'H2CO3': df,
    }
    df = pd.DataFrame()
    for label, df_p in dict_x.items():
        np_hist, bin_edges = np.histogram(df_p['roh0(ang)'], bins=200, density=True, range=(0.8,10))
        np_energy = plm.prob_to_deltag(np_hist, temperature=330)
        bin_center = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2
        df_tmp = pd.DataFrame({'roh0(ang)': bin_center, label:np_energy})
        x = df_tmp['roh0(ang)']
        y = df_tmp[label]
        #np_energy -= np.amin(np_energy)
        y -= min(y[(x>3.5) & (x<4.5)])
        df = pd.concat([df, y], axis=1)
    df = pd.concat([df, x], axis=1)

    df.to_csv('carbonic_roh_1d.csv', index=False)

def main():

    run()

main()
