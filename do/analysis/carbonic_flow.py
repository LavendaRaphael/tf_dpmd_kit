import pandas as pd
from tf_dpmd_kit import analysis
import json

def main():
    
    file_save = 'carbonic_flow.csv'

    list_file = [
        '../CC/carbonic/timelong.json',
        '../CT/carbonic/timelong.json',
        '../TT/carbonic/timelong.json',
    ]
    time_tot = 0
    for file in list_file:
        with open(file, 'r') as fp:
            time_tot += json.load(fp)['timelong(ps)']
    print('time_tot', time_tot)

    list_file = [
        '../CC/carbonic/volume.json',
        '../CT/carbonic/volume.json',
        '../TT/carbonic/volume.json',
    ]
    volume = 0
    for file in list_file:
        with open(file, 'r') as fp:
            volume += json.load(fp)['volume(ang3)']
    volume /= len(list_file)
    print('volume', volume)

    df = analysis.read_multidata(
        list_file = [
            '../CC/carbonic/carbonic_lifedata.csv',
            '../CT/carbonic/carbonic_lifedata.csv',
            '../TT/carbonic/carbonic_lifedata.csv',
        ]
    )
    df.dropna(inplace=True)
    print(df)

    list_flow = []
    col_0 = pd.Series({'event': 0.0})
    for idx, col in df.iterrows():
        if col_0['event'] == 1.0:
            list_flow.append([col_0['state'], col['state'], 1])
        elif col_0['event'] == 0.0:
            list_flow.append(['X', col['state'], 1])
        else:
            raise
        if col['event'] == 0.0:
            list_flow.append([col['state'], 'X', 1])
        col_0 = col
    df_flow = pd.DataFrame(list_flow, columns=['from', 'to', 'freq(ns-1)'])
    dfgb = df_flow.groupby(['from','to'])
    df_count = dfgb.count()/time_tot*1000

    Avogadro = 6.02214076e23
    molar = 1e27/Avogadro/volume
    ns2s = 1e9
    df_count['rate(M/s)'] = df_count['freq(ns-1)']*molar*ns2s
    print(file_save)
    print(df_count)
    df_count.to_csv(file_save)

main()
