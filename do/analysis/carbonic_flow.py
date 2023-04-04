import pandas as pd
from tf_dpmd_kit import analysis
import json

def main():
    
    list_timelong = [
        '../TT/carbonic/timelong.json',
        '../CT/carbonic/timelong.json',
        '../CC/carbonic/timelong.json',
    ]
    time_tot = 0
    for file in list_timelong:
        with open(file, 'r') as fp:
            dict_timelong = json.load(fp)
            time_tot += dict_timelong['timelong(ps)']

    df = analysis.data_merge(
        list_file = [
            '../TT/carbonic/carbonic_lifedata.csv',
            '../CT/carbonic/carbonic_lifedata.csv',
            '../CC/carbonic/carbonic_lifedata.csv',
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
            list_flow.append([col['state'], 'X', 1])
        else:
            raise
        if col['event'] == 0.0:
            list_flow.append(['X', col['state'], 1])
        col_0 = col
    df_flow = pd.DataFrame(list_flow, columns=['from', 'to', 'count'])
    dfgb = df_flow.groupby(['from','to'])
    ser_count = dfgb.count()/time_tot*1000
    print(ser_count)
main()
