import math
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def grid_plt(
    list2d_file: list[list],
    str_xlabel: str,
    str_ylabel: str,
    str_save: str=None,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    tup_colormap: tuple = None,
    bool_minzero: bool = False,
    bool_maxzero: bool = False,
    bool_minus: bool = False,
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    if tup_colormap:
        sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=tup_colormap[0], vmax=tup_colormap[1]))

    fig, ax = plt.subplots()
    for list_file in list2d_file:
        str_file = list_file[0]
        str_label = list_file[1]
        print(str_file)
        if not os.path.isfile(str_file):
            continue

        if tup_colormap:
            color = sm.to_rgba(int(str_label[:-1]))
        else:
            color = None

        np_data = np.loadtxt(str_file)
        np_data_y = np_data[:,1]
        if bool_minus:
            np_data_y *= -1

        if bool_minzero:
            np_data_y -= min(np_data_y)

        if bool_maxzero:
            np_data_y -= max(np_data_y)

        ax.plot( np_data[:,0], np_data_y, label=str_label, linewidth=2, color=color)
    
    ax.legend()
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel(str_ylabel)
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.set_size_inches(9, 7)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)

def colvar_plt(
    list_header: list,
    float_timescale: float = 1.0,
    tup_xlim: tuple = None,
    str_save: str = None,
    str_color: str = None,
    dict_label: dict = None,
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    int_nplot = len(list_header)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')
    if int_nplot==1:
        axs = [axs]

    with open('COLVAR', 'r') as file_open:
        str_line = file_open.readline()
        list_tmp = [str_line]
        list_field = str_line.split()[2:]
        int_count = 0
        for str_line in file_open:
            print(list_tmp)
            list_tmp.append(str_line)
            for str_line in file_open:
                if str_line[:9] != '#! FIELDS':
                    list_tmp.append(str_line)
                else:
                    break
            data = np.genfromtxt(list_tmp, dtype=None, names=list_field, invalid_raise=False)
            for int_i in range(int_nplot):
                str_header = list_header[int_i][0]
                if str_header not in data.dtype.names:
                    continue
                axs[int_i].scatter(data['time']*float_timescale, data[str_header], s=0.5, color=str_color)
            if str_line[0] != '#':
                break
            list_tmp = [str_line]
            list_field = str_line.split()[2:]
    for int_i in range(int_nplot):
        str_header = list_header[int_i][0]
        if str_header not in dict_label:
            str_label = str_header
        else:
            str_label = dict_label[str_header]
        axs[int_i].set_ylabel(str_label)
        axs[int_i].set_xlim(tup_xlim)
        if len(list_header[int_i]) > 1:
            axs[int_i].set_ylim(list_header[int_i][1])
    axs[-1].set_xlabel(('time (ns)'))
    if str_save:
        fig.set_size_inches(11, 5)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)

def hills_sum(
    str_min: str,
    str_max: str,
    str_bin: str,
    str_cv: str = None,
    str_outfile: str = None,
    str_in: str = None,
    str_log: str = None,
    str_hills: str = 'HILLS'
):

    float_T, float_KbT = get_temperature(str_in, str_log)

    if not str_outfile:
        str_outfile = str_cv+'_fes.'

    str_cmd = f'source {os.environ["HOME"]}/.config/.tianff &&'
    str_cmd += 'source ${homedir}/.local/bin/bashrc_plm.sh ;'
    str_cmd += f'plumed sum_hills --hills {str_hills} --stride 20000 --outfile {str_outfile} --min {str_min} --max {str_max} --bin {str_bin} --negbias --kt {float_KbT}'
    if str_cv:
        str_cmd += f' --idw {str_cv}'

    subprocess_results = subprocess.run( str_cmd, shell=True, check=True, text=True, executable='/bin/bash')

    print(subprocess_results.stdout)

def T2KbT(
    float_T: float,
) -> float:

    float_Avogadro = 6.02214076e23
    # J*K^-1
    float_Kb = 1.380649e-23
    # KJ*mol^-1
    float_KbT = float_Kb*float_Avogadro*float_T/1000.0
    return float_KbT

def get_temperature(
    str_in: str = '../plm.in',
    str_log: str = '../plm.log'
) -> float:

    if str_in:
        with open(str_in, 'r') as fp:
            for str_line in fp:
                if 'TEMP=' in str_line:
                    list_line = str_line.split()
                    break
            for str_key in list_line:
                if 'TEMP=' == str_key[:5]:
                    float_T = int(str_key[5:])
                    break
        if not float_T:
            raise
        float_T2KbT = T2KbT(float_T)

    if str_log:
        with open(str_log, 'r') as fp:
            for str_line in fp:
                if 'KbT' in str_line:
                    list_line = str_line.split()
                    float_KbT = float(list_line[2])
                    break
        if not float_KbT:
            raise

    if str_in and str_log:
        if float_T2KbT-float_KbT > 1e4:
            raise

    if not (str_in or str_log):
        raise

    return float_T, float_T2KbT

def deltag_to_pka(
    float_deltag: float,
    float_T: float = None,
    str_in: str = '../plm.in',
    str_log: str = '../plm.log'
) -> float:

    if float_T is None:
        float_T, float_KbT = get_temperature(str_plmin, str_plmlog)
    else:
        float_KbT = T2KbT(float_T)

    return float_deltag/(float_KbT*math.log(10))

def pka_to_deltag(
    float_pka: float,
    float_T: float = None,
    str_in: str = '../plm.in',
    str_log: str = '../plm.log'
) -> float:

    if float_T is None:
        float_T, float_KbT = get_temperature(str_plmin, str_plmlog)
    else:
        float_KbT = T2KbT(float_T)

    return float_pka*(float_KbT*math.log(10))

