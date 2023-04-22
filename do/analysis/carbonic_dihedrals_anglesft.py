import os
import pandas as pd
import numpy as np

def angle_sft(x):

    if x < -np.pi/2:
        return x + 2*np.pi
    else:
        return x

def run():

    file_data = 'carbonic.product.csv'
    file_save = 'carbonic_dihedrals.csv'

    print(file_data)
    df = pd.read_csv(file_data, index_col=0)

    df = df[df['dihedral0(rad)'].notnull()]
    df = df.loc[:, ['dihedral0(rad)', 'dihedral1(rad)']]
    df = df.applymap(angle_sft)

    print(file_save)
    print(df)
    df.to_csv(file_save)

