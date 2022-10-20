import json
import dpdata

def check_locate(
    int_id: int,
    str_data: str
):
    with open(str_data, 'r') as fp:
        dict_json = json.load(fp)
    dict_system = dict_json['dict_system']

    int_tmp = int_id
    for str_system in dict_system:
        int_num = dict_system[str_system]
        if int_tmp < int_num:
            return str_system, int_tmp
        else:
            int_tmp -= int_num
    print(int_id)
    raise

def from_input():
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

def from_dptest():
    with open('dptest.e.out', 'r') as fp:
        str_system = fp.readline().split()[1]
        for str_line in fp:
            int_tmp = 1
            for str_line in fp:
                if str_line[0] == '#':
                    print(str_system, int_tmp)
                    str_system = str_line.split()[1]
                    break
                else:
                    int_tmp += 1
