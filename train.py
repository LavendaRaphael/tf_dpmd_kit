import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import dpdata
from tf_dpmd_kit import plot

def data_check_locate(
    int_id: int,
    str_json_data: str
):
    with open(str_json_data, 'r') as fp:
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

def datastatus_from_input():
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

def datastatus_from_dptest():
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

def dptest_plt(
    str_file: str,
    int_natoms: int = None,
    str_plt_type: str = 'dft_dnn',
    tup_hist_xlim: tuple = None,
    float_lw: float = None,
) -> None:

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

    fig, ax = plt.subplots()

    float_rmse_e, float_rmse_f, float_rmse_v = get_rmse()

    if (str_mode=='e'):
        # per atom
        np_data_new = np_data / int_natoms
        # eV to meV
        float_rmse_e *= 1000
        np_data_new *= 1000
        np_data_new -= np.average( np_data_new[:,0] )

        tup_dft_dnn_label = ('DFT energy (meV/atom)', 'DNN energy (meV/atom)')
        str_label = f'Energy RMSE = {float_rmse_e:.3f} meV/atom'
        str_title = 'Energy'
        str_inset_xlabel = 'error'
        str_hist_xlabel = r'E$_{DP}$-E$_{DFT}$ (meV/atom)'

    elif (str_mode=='f'):
        np_data_new = np_data.reshape((np_data.shape[0]*3, 2), order='F')

        tup_dft_dnn_label = ('DFT force (eV/Å)', 'DNN force (eV/Å)')
        str_label = f'Force RMSE = {float_rmse_f:.3f} eV/Å'
        str_title = 'Force'
        str_inset_xlabel = 'error'
        str_hist_xlabel = r'F$_{DP}$-F$_{DFT}$ (eV/Å)'

    del_data = np_data_new[:,1] - np_data_new[:,0]
    if (str_plt_type == 'dft_dnn' or str_plt_type == 'inset'):
        ax.scatter(
            np_data_new[:,0],
            np_data_new[:,1],
            #label =str_label,
            edgecolors='none', 
            s=1.5,
        )
        ax.set_xlabel(tup_dft_dnn_label[0])
        ax.set_ylabel(tup_dft_dnn_label[1])
        ax.axline([0, 0], [1, 1], color='black', linestyle='--', lw=float_lw)
        if str_plt_type == 'inset':
            axin = ax.inset_axes((0.7,0.2,0.28,0.22))
            axin.hist(
                del_data,
                bins = 'auto',
                density = True
            )
            axin.set_xlim(tup_hist_xlim)
            axin.set_xlabel(str_inset_xlabel)
            axin.set_ylabel('Prob. Dens.')

    elif (str_plt_type == 'hist'):
        ax.hist(
            del_data,
            #label = str_label,
            bins = 'auto',
            density = True
        )
        ax.set_xlabel(str_hist_xlabel)
        ax.set_xlim(tup_hist_xlim)
        ax.set_ylabel('Probability Density')
    else:
        raise

    plot.set_lw(ax, float_lw)

    ax.legend(
        title = str_title,
        frameon = False,
        loc = 'upper left'
    )

    return fig

def get_rmse():
    with open('log', 'r') as fp:
        for str_line in fp:
            if 'Energy RMSE/Natoms' in str_line:
                float_rmse_e = float(str_line.split()[-2])
            elif 'Force  RMSE' in str_line:
                float_rmse_f = float(str_line.split()[-2])
            elif 'Virial RMSE/Natoms' in str_line:
                float_rmse_v = float(str_line.split()[-2])
    return float_rmse_e, float_rmse_f, float_rmse_v

def def_plt_lcurve(
    str_save: str = None,
    tup_size: tuple = None,
) -> None:

    mpl.rcParams['font.size']=15
    mpl.rcParams['font.family']='sans-serif'
    mpl.rcParams['font.sans-serif']=["Arial"]

    fig, ax = plt.subplots()

    data = np.genfromtxt("lcurve.out", names=True)
    ax.scatter(data['step'], data['rmse_e_trn'], label='Energy', s=0.5)
    ax.scatter(data['step'], data['rmse_f_trn'], label='Force', s=0.5)
    ax.legend(loc='lower left')
    ax.set_xlabel('Step')
    ax.set_ylabel('Loss')
    ax.set_yscale('log')
    ax.set_xlim((None,None))
    ax.grid()
