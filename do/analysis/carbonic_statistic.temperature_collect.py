import pandas as pd

def run():

    dict_data = {
        290: '../290K/carbonic/carbonic_statistic.csv',
        310: '../310K/carbonic/carbonic_statistic.csv',
        330: '../330K/carbonic/carbonic_statistic.csv',
        350: '../350K/carbonic/carbonic_statistic.csv',
    }

    list_header = pd.read_csv(list(dict_data.values())[0], index_col=0).columns
    print(list_header)

    dict_header = {
        'timemean(ps)': 'timemean',
        'count(ns-1)': 'count'
    }

    for header in list_header:
        print(header)

        df_save = pd.DataFrame()
        for temperature, file_data in dict_data.items():
            df_data = pd.read_csv(file_data, index_col=0)
            df_x = df_data[header].rename(temperature).to_frame().T
            df_save = pd.concat([df_save, df_x])

        df_save.index.name = 'temperature(K)'
        if header in dict_header:
            name = dict_header[header]
        else:
            name = header
        file_save = f'carbonic_statistic.temperature_{name}.csv'
        print(file_save)
        print(df_save)
        df_save.to_csv(file_save)

run()
