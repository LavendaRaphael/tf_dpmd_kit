import math
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from tf_dpmd_kit import plot
import pandas as pd

def colvar_hist_plt(
    str_header: str,
    str_label: str,
    str_save: str = None,
    tup_timerange: tuple = None,
    tup_xlim: tuple = None,
):

    fig, ax = plt.subplots()

    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    np_data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
    #print(data.dtype)

    if tup_timerange is None:
        np_data_new = data[str_header]
    else:
        np_data_time = np_data['time']
        list_tof = (np_data_time >= tup_timerange[0]) & (np_data_time <= tup_timerange[1])
        np_data_time_new = np_data_time[list_tof]
        print(np_data_time_new[0],np_data_time_new[-1])
        np_data_new = np_data[str_header][ list_tof ]

    ax.hist(
        np_data_new,
        bins = 'auto',
        density = True
    )

    float_std = np.std(np_data_new)
    float_mean = np.mean(np_data_new)

    ax.legend(
        title = f'MEAN = {float_mean:.3f}\nSTD = {float_std:.3f}'
    )
    ax.set_xlim(tup_xlim)
    ax.set_xlabel(str_label)
    ax.set_ylabel('Probability Density')

    return fig, ax

def fes_integral(
    str_file: str,
    tup_xrange: tuple = None,
    float_T: float = None,
    str_in: str = None,
    str_log: str = None
) -> float:

    np_data = np.loadtxt(str_file)    
    np_data[:,1] -= min(np_data[:,1])

    if not float_T:
        float_T, float_KbT = get_temperature( str_in, str_log )
    else:
        float_KbT = T2KbT(float_T)
    
    np_data[:,1] = np.exp(-np_data[:,1]/float_KbT)

    if tup_xrange is None:
        np_data_trapz = np_data
    else:
        float_xlow = tup_xrange[0]
        float_xup = tup_xrange[1]
        if float_xup == '1M':
            float_Avogadro = 6.02214076e23
            float_xup = 1e27/float_Avogadro
            float_xup = (float_xup/(4/3*math.pi))**(1.0/3.0)

        np_indices = np.searchsorted(np_data[:,0], [float_xlow, float_xup])
        slice_range = slice(np_indices[0], np_indices[1])
        np_data_trapz = np_data[slice_range, :]

    float_integral = np.trapz(
        np_data_trapz[:, 1],
        x = np_data_trapz[:, 0]
    )

    float_integral = -np.log(float_integral) * float_KbT
    return float_integral

def fes_1M_correct(
    float_volume: float,
    float_T: float
) -> float:

    float_Avogadro = 6.02214076e23
    float_volume_1M = 1e27/float_Avogadro
    float_KbT = T2KbT(float_T)
    #print(float_volume_1M/float_volume,'mol/L')
    float_correct = -float_KbT*math.log(float_volume_1M/float_volume)

    return float_correct

def get_pka(
    str_file: str,
    dict_coef: dict,
    tup_srange_tot: tuple,
    float_volume: float = None,
    float_T: float = None,
    str_in: str = None,
    str_log: str = None
) -> (float, float):

    '''Compute pka.

    Args:
        dict_coef: 
            {
                (s0, s1): c0,
                (s1, s2): c1,
                ...
            }
        float_volume: volume
            V
        tup_srange_tot:
            (sl, sr)
        float_T: temperature
            T

    Returns:
          DeltaF 
        = 
          c0*[-kB*T*ln(P(s0,s1)/P(sl,sr))]
        + c1*[-kB*T*ln(P(s1,s2)/P(sl,sr))]
        + [-kB*T*ln(V_x/V)]

    '''

    if not float_T:
        float_T, float_KbT = get_temperature( str_in, str_log )

    float_h_tot = fes_integral(
        str_file = str_file,
        tup_xrange = tup_srange_tot,
        float_T = float_T,
    )
 
    float_deltag = 0
    for tup_srange, float_coef in dict_coef.items():
        float_h = fes_integral(
            str_file = str_file,
            tup_xrange = tup_srange,
            float_T = float_T,
        )
        float_deltag += float_coef*(float_h-float_h_tot)

    if not (float_volume is None):
        float_correct = fes_1M_correct(
            float_volume = float_volume,
            float_T = float_T
        )
        float_deltag += float_correct

    float_pka = deltag_to_pka(
        float_deltag = float_deltag,
        float_T = float_T,
    )

    return float_deltag, float_pka

def get_pka_time(
    dict_file: dict,
    dict_coef: dict,
    float_volume: float = None,
    tup_srange_tot: tuple = None,
    str_save: str = None,
    float_T: float = None,
    str_in: str = None,
    str_log: str = None
):

    int_nfile = len(dict_file)
    np_pka = np.zeros(shape=(int_nfile), dtype=[('time', 'f4'), ('pka', 'f4')])
    np_deltag = np.zeros(shape=(int_nfile), dtype=[('time', 'f4'), ('deltag', 'f4')])
    if not float_T:
        float_T, float_KbT = get_temperature( str_in, str_log )
    print(float_T)
    for int_i,int_key in enumerate(dict_file):
        if not os.path.isfile(dict_file[int_key]):
            continue
        float_deltag, float_pka = get_pka(
            str_file = dict_file[int_key],
            dict_coef = dict_coef,
            tup_srange_tot = tup_srange_tot,
            float_T = float_T,
            float_volume = float_volume
        )
        np_pka[int_i] = (int_key, float_pka)
        np_deltag[int_i] = (int_key, float_deltag)

    if str_save:
        str_save_pka = f'{str_save}.pka.csv'
        print(str_save_pka)
        np.savetxt(str_save_pka, np_pka, header=' '.join(np_pka.dtype.names))
        str_save_deltag = f'{str_save}.deltag.csv'
        print(str_save_deltag)
        np.savetxt(str_save_deltag, np_deltag, header=' '.join(np_deltag.dtype.names))

def colvar_reader(
    list_data: list,
):

    df_data = pd.DataFrame()
    for str_file in list_data:
        with open(str_file, 'r') as file_open:
            str_line = file_open.readline()
            list_field = str_line.split()[2:]
        df_tmp = pd.read_csv(str_file, sep=' ', comment='#', names=list_field)
        df_data = pd.concat([df_data, df_tmp])

    return df_data

def colvar_plt(
    dict_header: dict,
    list_data: list = None,
    str_xlabel: str = 'Time (ps)',
    float_timescale: float = 1.0,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    str_color: str = 'black',
    dict_color: dict = None,
    str_title: str = None,
    legend_loc: str = None,
    dict_label: dict = None,
    bool_scatter: bool = True,
    float_lw: float = None,
) -> None:

    if list_data is None:
        list_data = ['COLVAR']

    int_nplot = len(dict_header)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')

    if int_nplot==1:
        axs = [axs]

    df_data = colvar_reader(list_data)

    for int_i,str_header in enumerate(dict_header):
        if bool_scatter:
            axs[int_i].scatter(df_data['time']*float_timescale, df_data[str_header], s=1, color=str_color, edgecolors='none')
            if not(dict_color is None):
                for tup_range, str_color_i in dict_color.items():
                    df_tmp = df_data.iloc[tup_range[0]:tup_range[1]]
                    if not(tup_range in dict_label):
                        dict_label[tup_range] = None
                    axs[int_i].scatter(df_tmp['time']*float_timescale, df_tmp[str_header], s=1, color=str_color_i, edgecolors='none', label=dict_label[tup_range])
        else:
            axs[int_i].plot(df_data['time']*float_timescale, df_data[str_header], color=str_color, linewidth=float_lw)

        str_ylabel = dict_header[str_header]
        axs[int_i].set_ylabel(str_ylabel)

        plot.set_lw( axs[int_i], float_lw )

        axs[int_i].legend(
            loc = legend_loc,
            frameon = False,
            markerscale = 3,
            labelspacing = 0.3,
        )
    axs[0].set_xlim(tup_xlim)
    axs[0].set_ylim(tup_ylim)
    if str_xlabel:
        axs[-1].set_xlabel(str_xlabel)

    return fig, axs

def colvar_plt_old(
    dict_header: dict,
    list_data: list = None,
    str_xlabel: str = 'Time (ps)',
    float_timescale: float = 1.0,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    str_color: str = None,
    str_title: str = None,
    str_legeng_loc: str = None,
    bool_scatter: bool = True,
    float_lw: float = None,
) -> None:

    if list_data is None:
        list_data = ['COLVAR']

    int_nplot = len(dict_header)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')

    if int_nplot==1:
        axs = [axs]

    for str_file in list_data:
        with open(str_file, 'r') as file_open:
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
                for int_i,str_header in enumerate(dict_header):
                    if str_header not in data.dtype.names:
                        continue
                    #for tup_
                    if bool_scatter:
                        axs[int_i].scatter(data['time']*float_timescale, data[str_header], s=1, color=str_color, edgecolors='none')
                    else:
                        axs[int_i].plot(data['time']*float_timescale, data[str_header], color=str_color, linewidth=float_lw)

                if str_line[0] != '#':
                    break
                list_tmp = [str_line]
                list_field = str_line.split()[2:]

    for int_i,str_header in enumerate(dict_header):
        str_ylabel = dict_header[str_header]
        axs[int_i].set_ylabel(str_ylabel)
        plot.set_lw( axs[int_i], float_lw )

    if str_title:
        axs[0].legend(
            loc = str_legeng_loc,
            title = str_title,
            frameon = False,
        )
    axs[0].set_xlim(tup_xlim)
    axs[0].set_ylim(tup_ylim)
    if str_xlabel:
        axs[-1].set_xlabel(str_xlabel)

    return fig, axs

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

