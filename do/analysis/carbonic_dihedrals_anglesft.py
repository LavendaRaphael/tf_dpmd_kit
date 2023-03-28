import pandas as pd
import numpy as np

def angle_sft(x):

    if x < -np.pi/2:
        return x + 2*np.pi
    else:
        return x

def run(
    list_file: list,
    file_save: str,
):

    df_new = None
    for str_file in list_file:
        print(str_file)
        df_tmp = pd.read_csv(str_file)
        if df_new is None:
            df_new = df_tmp
        else:
            df_new = pd.concat([df_new, df_tmp], ignore_index=True)
    print(df_new)

    df_new = df_new[df_new['dihedral0(rad)'].notnull()]
    df_new = df_new.loc[:, ['dihedral0(rad)', 'dihedral1(rad)']]
    df_new = df_new.applymap(angle_sft)
    print(df_new)

    df_new.to_csv(file_save, index=False)

run(
    list_file = [
        '../TT/carbonic/carbonic.product.csv',
        '../CT/carbonic/carbonic.product.csv',
        '../CC/carbonic/carbonic.product.csv',
    ],
    file_save = 'carbonic_dihedrals.csv'
)
