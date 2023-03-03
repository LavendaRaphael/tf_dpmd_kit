import json
import numpy as np

def run(
    dict_file: dict,
    str_save: str,
) -> None:
    
    np_data = np.zeros( shape=(len(dict_file),4) )
    for int_id, int_key in enumerate(dict_file):
        with open(dict_file[int_key], 'r') as fp:
            dict_tmp = json.load(fp)
        np_data[int_id, 0] = int_key
        np_data[int_id, 1] = dict_tmp['f_mean']
        np_data[int_id, 2] = dict_tmp['f_sem']
        np_data[int_id, 3] = dict_tmp['f_std']

    np.savetxt(
        fname = str_save,
        X = np_data,
        header = 'temperature f_mean f_sem f_std'
    )

list_tmp = [
    'o_0.o_w.OhWaterAngle', 
    'o_0.o_w.OhWaterDistance',
    'o_1.o_w.OhWaterAngle', 
    'o_1.o_w.OhWaterDistance',
]

for str_tmp in list_tmp:

    str_file = f'hbondslength.{str_tmp}.ave.json'
    run(
        dict_file = {
            280: f'../280K/hbondslength/{str_file}',
            290: f'../290K/hbondslength/{str_file}',
            300: f'../300K/hbondslength/{str_file}',
            310: f'../310K/hbondslength/{str_file}',
            320: f'../320K/hbondslength/{str_file}',
            330: f'../330K/hbondslength/{str_file}',
        },
        str_save = f'hbondslength.{str_tmp}.ave.csv',
    )
    
