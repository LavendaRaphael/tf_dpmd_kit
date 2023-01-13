import os
import numpy as np

def hbonds_collect(
    dict_file: dict,
    str_save: str,
) -> None:
    
    with open(list(dict_file.values())[0], 'r') as fp:
        list_header = fp.readline().split()[1:]

    for str_header in list_header:
        np_data = np.zeros( shape=(len(dict_file),3) )
        for int_id, int_key in enumerate(dict_file):
            np_tmp = np.genfromtxt(dict_file[int_key], names=True)
            np_data[int_id, 0] = int_key
            np_data[int_id, 1] = np_tmp[str_header][0]
            np_data[int_id, 2] = np_tmp[str_header][2]

        np.savetxt(
            fname = f'{str_save}.{str_header}.csv',
            X = np_data,
            header = f'temperature {str_header} error'
        )

list_pair = [
    'h_0_1.o_w',
]

for str_pair in list_pair:
    dict_file = {
        280: f'280K/hbonds_new/hbonds.{str_pair}.ave.csv',
        290: f'290K/hbonds_new/hbonds.{str_pair}.ave.csv',
        300: f'300K/hbonds_new/hbonds.{str_pair}.ave.csv',
        310: f'310K/hbonds_new/hbonds.{str_pair}.ave.csv',
        320: f'320K/hbonds_new/hbonds.{str_pair}.ave.csv',
        330: f'330K/hbonds_new/hbonds.{str_pair}.ave.csv',
    }
    
    hbonds_collect(
        dict_file = dict_file,
        str_save = f'hbonds.{str_pair}.ave.temperature.new',
    )
    
