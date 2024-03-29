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

def dptest_parity_plt(
    ax,
    file_out: str,
    int_natoms: int = None,
    float_lw: float = None,
    list_ticks: list = None,
    file_log: str = 'log',
    bool_rmse: bool = False,
    energy_sep: bool = False,
) -> None:

    with open(file_out, 'r') as file_open:
        list_line = file_open.readline().split()
        if list_line[-1] == 'pred_e':
            mode = 'e'
        elif list_line[-1] == 'pred_fz':
            mode = 'f'
        else:
            print(list_line)
            raise

    np_data = np.loadtxt(file_out)

    if (mode=='e'):
        # per atom
        np_data_new = np_data / int_natoms
        # eV to meV
        np_data_new *= 1000
        np_data_new -= np.average( np_data_new[:,0] )
        if energy_sep:
            # two data set, the energy is sperate
            list_tof = (np_data_new[:,0] > 0)
            np_data_0 = np_data_new[list_tof]
            np_data_0 -= np.average(np_data_0)
            np_data_1 = np_data_new[~list_tof]
            np_data_1 -= np.average(np_data_1)
            np_data_new = np.concatenate((np_data_0, np_data_1), axis=0)
        
        xlabel = 'DFT energy (meV/atom)'
        ylabel = 'DP energy (meV/atom)'
        title = 'Energy'

    elif (mode=='f'):
        np_data_new = np_data.reshape((np_data.shape[0]*3, 2), order='F')
        unit = 'eV/Å'
        xlabel = 'DFT force (eV/Å)'
        ylabel = 'DP force (eV/Å)'
        title = 'Force'

    ax.axline([0, 0], [1, 1], linestyle='--', lw=float_lw, color='tab:orange')
    ax.scatter(
        np_data_new[:,0],
        np_data_new[:,1],
        edgecolors='none', 
        s=1.5,
        rasterized=True,
        color = 'tab:blue',
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    xmax = np.max(np.absolute(np_data_new))*1.1
    tup_lim = (-xmax, xmax)
    ax.set_ylim(tup_lim)
    ax.set_xlim(tup_lim)
    #ax.set_aspect(1)
    if not (list_ticks is None):
        ax.set_xticks(list_ticks)
        ax.set_yticks(list_ticks)

    plot.set_lw(ax, float_lw)

    ax.legend(
        title = title,
        frameon = False,
        loc = 'upper left'
    )
    if bool_rmse:
        rmse_e, rmse_f, rmse_v = get_rmse(file_log)
        if mode == 'e':
            rmse_e *= 1000
            s = f'RMSE = {rmse_e:.2f} meV/atom'
        elif mode == 'f':
            s = f'RMSE = {rmse_f:.3f} eV/Å'
        ax.text(
            x = 0.95,
            y = 0.05,
            s = s,
            ha = 'right',
            va = 'bottom',
            transform = ax.transAxes,
        )

def dptest_hist_plt(
    ax,
    file_out: str,
    int_natoms: int = None,
    xlim: tuple = None,
    bins: int = 'auto',
    float_lw: float = None,
) -> None:

    with open(file_out, 'r') as file_open:
        list_line = file_open.readline().split()
        if list_line[-1] == 'pred_e':
            mode = 'e'
        elif list_line[-1] == 'pred_fz':
            mode = 'f'
        else:
            print(list_line)
            raise

    np_data = np.loadtxt(file_out)

    if (mode=='e'):
        # per atom
        np_data_new = np_data / int_natoms
        # eV to meV
        np_data_new *= 1000
        np_data_new -= np.average( np_data_new[:,0] )
        xlabel = r'E$_{DP}$-E$_{DFT}$ (meV/atom)'
        title = 'Energy'

    elif (mode=='f'):
        np_data_new = np_data.reshape((np_data.shape[0]*3, 2), order='F')
        xlabel = r'F$_{DP}$-F$_{DFT}$ (eV/Å)'
        title = 'Force'

    del_data = np_data_new[:,1] - np_data_new[:,0]
    ax.hist(
        del_data,
        bins = 'auto',
        density = True
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Probability Density')
    plot.set_lw(ax, float_lw)

    ax.legend(
        title = title,
        frameon = False,
        loc = 'upper left'
    )

def get_rmse(
    file_log: str = 'log',
):
    with open(file_log, 'r') as fp:
        for str_line in fp:
            if 'Energy RMSE/Natoms' in str_line:
                float_rmse_e = float(str_line.split()[-2])
            elif 'Force  RMSE' in str_line:
                float_rmse_f = float(str_line.split()[-2])
            elif 'Virial RMSE/Natoms' in str_line:
                float_rmse_v = float(str_line.split()[-2])
    return float_rmse_e, float_rmse_f, float_rmse_v

def plt_lcurve(ax):

    data = np.genfromtxt("lcurve.out", names=True)
    ax.scatter(data['step'], data['rmse_e_trn'], label='Energy (eV)', s=0.5)
    ax.scatter(data['step'], data['rmse_f_trn'], label='Force (eV/Å)', s=0.5)
    ax.legend(loc='lower left')
    ax.set_xlabel('Step')
    ax.set_ylabel('Error')
    ax.set_yscale('log')
    ax.set_xlim((None,None))
    ax.grid()
