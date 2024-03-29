import pandas as pd

def run():

    dict_dir = {
        290: '../290K/carbonic/',
        310: '../310K/carbonic/',
        330: '../330K/carbonic/',
        350: '../350K/carbonic/',
    }

    list_key = []
    list_df = []
    for temperature, dir_x in dict_dir.items():
        df_data = pd.read_csv(dir_x+'carbonic_statistic.csv', index_col=0)
        df_life = pd.read_csv(dir_x+'carbonic_lifetime.csv', index_col=0)
        df_data = pd.concat([df_data, df_life], axis=1)

        list_key.append(temperature)
        list_df.append(df_data)
    df_save = pd.concat(list_df, keys=list_key, names=['temperature(K)'])
    df_save = df_save[['frac', 'frac_sem', 'freq(ns-1)', 'freq(ns-1)_sem', 'freqfrac', 'freqfrac_sem', 'lifetime(ps)', 'lower', 'upper', 'rate(M/s)', 'rate(M/s)_sem']]
    file_save = f'carbonic_statistic.temperature.csv'
    print(file_save)
    print(df_save)
    df_save.to_csv(file_save)

run()
