import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import json

def check_locate(
    int_id: int,
    str_data: str
):
    with open(str_data, 'r') as fp:
        dict_json = json.load(fp)
    dict_system = dict_json['dict_system']

    for str_system in dict_system:
        int_num = dict_system[str_system]
        if int_id < int_num:
            return str_system, int_id
        else:
            int_id -= int_num
    raise

def def_plt(
    float_rmse: float,
    str_file: str,
    int_natoms: int,
    str_save: str = None,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    int_bins: int = 'auto',
    str_data: str = None
) -> None:
    
    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    fig, ax = plt.subplots()

    with open(str_file, 'r') as file_open:
        list_line = file_open.readline().split()
        if list_line[-1] == 'pred_e':
            str_mode = 'e'
        elif list_line[-1] == 'pred_fz':
            str_mode = 'f'
        else:
            print(list_line)
            raise
    np_data = np.loadtxt(str_file)
    if (str_mode=='e'):
        # per atom
        np_data /= int_natoms
        # eV to meV
        np_data *= 1000
        del_data = np_data[:,1] - np_data[:,0]
        ax.set_xlabel( r'E$_{DP}$-E$_{DFT}$ (meV/atom)')
        str_label = f'Energy RMSE = {float_rmse:.3f} meV/atom'
        if str_data:
            for int_i in range(len(del_data)):
                if abs(del_data[int_i]) > 1.5:
                    str_locate, int_locate = check_locate(int_i, str_data)
                    print(f'{str_locate} {int_locate} {del_data[int_i]:.1f}')

    elif (str_mode=='f'):
        del_data_xyz = np_data[:,3:6] - np_data[:,0:3]
        del_data = np.linalg.norm( del_data_xyz, axis=1 )
        del_data *= 1000
        ax.set_xlabel(r'|F$_{DP}$-F$_{DFT}$| (meV/Å)')
        str_label = f'Force RMSE = {float_rmse:.1f} meV/Å'
        if str_data:
            for int_i in range(len(del_data)):
                if del_data[int_i] > 1000:
                    int_id = int_i // int_natoms
                    str_locate, int_locate = check_locate(int_id, str_data)
                    print(f'{str_locate} {int_locate} {del_data[int_i]:.1f}')

    #float_std = np.std(del_data)
    #float_mean = np.mean(del_data) 
    ax.hist(
        del_data,
        label = str_label,
        bins = int_bins,
        #density = True
    )

    ax.legend()
    ax.set_ylabel('Probabilty Density')
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

#'''
def_plt(
    str_file = 'dptest.e.out',
    int_natoms = 384,
    float_rmse = 0.412,
    #tup_xlim = (-1.5,1.5),
    #tup_ylim = (0,10),
    str_save = 'dptest.e.pdf',
    str_data = 'data.json'
)
#'''
#'''
def_plt(
    str_file = 'dptest.f.out',
    float_rmse = 59.9,
    int_natoms = 384,
    tup_xlim = (0,250),
    #tup_ylim = (0,10),
    str_save = 'dptest.f.pdf',
    str_data = 'data.json'
)
#'''

plt.show()
