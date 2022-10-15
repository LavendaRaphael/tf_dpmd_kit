import json
import dpdata

with open('input.json', 'r') as fp:
    dict_dpin = json.load(fp)
list_system = dict_dpin['training']['training_data']['systems']
dict_system = {}
int_sum = 0
for str_system in list_system:
    dp_sys = dpdata.LabeledSystem(file_name=str_system, fmt='deepmd/npy')
    int_num = len(dp_sys)
    dict_system[str_system] = int_num
    int_sum += int_num
    print(str_system, int_num)

dict_json = {
    'int_sum': int_sum,
    'dict_system': dict_system
}
with open('data.json', 'w') as fp:
    json.dump( dict_json, fp, indent=4 )
