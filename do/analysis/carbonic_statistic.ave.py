import pandas as pd
import os
import numpy as np

def main():

    list_file = [
        '../TT/carbonic/carbonic_statistic.csv',
        '../CT/carbonic/carbonic_statistic.csv',
        '../CC/carbonic/carbonic_statistic.csv',
    ]


    list_df = []
    list_key = []
    for file in list_file:
        df_tmp = pd.read_csv(file, index_col=0)
        list_df.append(df_tmp)
        list_key.append(file)
    df_data = pd.concat(list_df, keys=list_key)
    df_data = df_data.fillna(0)
    print(df_data)
    df_mean = df_data.groupby(level=1).mean()
    print(df_mean)
    df_std = df_data.groupby(level=1).std()
    print(df_std)
main()
