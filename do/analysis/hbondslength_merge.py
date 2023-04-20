import pandas as pd
import matplotlib.pyplot as plt

def run(
    dict_data: dict,
    dict_range: dict,
    str_header: str,
    str_save: str,
):

    df_data = pd.DataFrame()

    for key, str_data in dict_data.items():
        print(str_data)
        tmp_data = pd.read_csv(
            str_data,
            sep = ' ')
        tup_range = dict_range[key]
        tmp_data = tmp_data[tup_range[0]: tup_range[1]]
        print(tmp_data)
        df_data = pd.concat([df_data, tmp_data])
    print(str_save)
    print(df_data)
    df_data[str_header].to_csv(str_save, index=False)

list_header = ['OhWaterAngle', 'OhWaterDistance']

for str_header in list_header:
    run(
        dict_data = {
            '1': '../../../H2CO3_CC_H2O_126/280K/hbondslength/hbondslength.o_0.o_w.csv',
            '2': '../../../H2CO3_CC_H2O_126/280K/hbondslength/hbondslength.o_2.o_w.csv',
            '3': 'hbondslength.o_0.o_w.csv',
            '4': 'hbondslength.o_2.o_w.csv',
        },
        dict_range = {
            '1': (100000, 300000),
            '2': (100000, 300000),
            '3': (100000, 400000),
            '4': (100000, 400000),
        },
        str_header = str_header,
        str_save = f'hbondslength.o_0_2.o_w.{str_header}.csv',
    )
