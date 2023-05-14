import os
import pandas as pd
import numpy as np

def product(df):

    return df[100000:]

def roh_order(df):

    bool_x = df['roh0(ang)'] < df['roh1(ang)']
    df_0 = df.loc[bool_x]
    df_0.rename(columns={'dihedral0(rad)': 'dihedral1(rad)', 'dihedral1(rad)': 'dihedral0(rad)', 'roh0(ang)': 'roh1(ang)', 'roh1(ang)': 'roh0(ang)'}, inplace=True, errors='raise')
    return pd.concat([df[~bool_x], df_0])

def dh_sft(df):

    bool_x = df['dihedral0(rad)'] < -np.pi/2
    df.loc[bool_x, 'dihedral0(rad)'] += 2*np.pi

    bool_x = df['dihedral1(rad)'] < -np.pi/2
    df.loc[bool_x, 'dihedral1(rad)'] += 2*np.pi

    return df

def run():

    file_data = 'carbonic.csv'
    file_save = 'carbonic.product.csv'

    df = pd.read_csv(file_data, index_col=0)
    print(file_data)
    print(df)

    df = product(df)

    df = roh_order(df)

    df = dh_sft(df)

    print(file_save)
    print(df)
    df.to_csv(file_save)

def main():

    dir_list = [
        '290K/CC',
        '290K/CT',
        '290K/TT',
        '310K/CC',
        '310K/CT',
        '310K/TT',
        '330K/CC',
        '330K/CT',
        '330K/TT',
        '350K/CC',
        '350K/CT',
        '350K/TT',
    ]
    dir_cwd = os.getcwd()
    for dir_i in dir_list:
        dir_i = os.path.join('..', dir_i, 'carbonic')
        print(dir_i)
        os.chdir(dir_i)

        run()

        os.chdir(dir_cwd)

main()
