import pandas as pd
import numpy as np

def main():

    list_file = [
        '../TT/carbonic/carbonic_statistic.csv',
        '../CT/carbonic/carbonic_statistic.csv',
        '../CC/carbonic/carbonic_statistic.csv',
    ]

    file_mean = 'carbonic_statistic.mean.csv'
    file_sem = 'carbonic_statistic.sem.csv'

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
    print(file_mean)
    print(df_mean)
    df_mean.to_csv(file_mean)

    df_sem = df_data.groupby(level=1).sem()
    print(df_sem)
    print(file_sem)
    df_sem.to_csv(file_sem)
    
main()
