import pandas as pd
import json

def run(
    str_data,
    str_save,
):

    df_data = pd.read_csv(str_data)
    dict_save = {}
    dict_save['f_mean'] = df_data.mean()[0]
    dict_save['f_std'] = df_data.std()[0]
    dict_save['f_sem'] = df_data.sem()[0]
    print(dict_save)
    with open(str_save, 'w') as fp:
        json.dump(dict_save, fp, indent=4)

run(
    str_data = 'hbondslength.o_0_2.o_w.OhWaterDistance.csv',
    str_save = 'hbondslength.o_0_2.o_w.OhWaterDistance.ave.json',
)
run(
    str_data = 'hbondslength.o_0_2.o_w.OhWaterAngle.csv',
    str_save = 'hbondslength.o_0_2.o_w.OhWaterAngle.ave.json',
)
