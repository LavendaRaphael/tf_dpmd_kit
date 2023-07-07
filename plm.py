import math
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from tf_dpmd_kit import plot
import pandas as pd

def colvar_hist_plt(
    ax,
    header: str,
    file_data: str = 'COLVAR',
    timerange: tuple = None,
):

    with open(file_data, 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    np_data = np.genfromtxt(file_data, dtype=None, names=list_header)

    if timerange is None:
        np_data_new = data[header]
    else:
        np_data_time = np_data['time']
        list_tof = (np_data_time >= timerange[0]) & (np_data_time <= timerange[1])
        np_data_time_new = np_data_time[list_tof]
        np_data_new = np_data[header][ list_tof ]

    ax.hist(
        np_data_new,
        bins = 'auto',
        density = True
    )

    std = np.std(np_data_new)
    mean = np.mean(np_data_new)

    ax.legend(
        title = f'MEAN = {mean:.3f}\nSTD = {std:.3f}',
        frameon = False,
    )
    ax.set_ylabel('Probability Density')

def fes_integral(
    data_x,
    data_y,
    tup_xrange: tuple = None,
    temperature: float = None,
) -> float:

    data_y -= min(data_y)

    kBT = T2kBT(temperature)
    
    data_y = np.exp(-data_y/kBT)

    if tup_xrange is None:
        data_trapz_x = data_x
        data_trapz_y = data_y
    else:
        xlow = tup_xrange[0]
        xup = tup_xrange[1]
        if xup == '1M':
            Avogadro = 6.02214076e23
            xup = 1e27/Avogadro
            xup = (xup/(4/3*math.pi))**(1.0/3.0)

        np_indices = np.searchsorted(data_x, [xlow, xup])
        slice_range = slice(np_indices[0], np_indices[1])
        data_trapz_x = data_x[slice_range]
        data_trapz_y = data_y[slice_range]

    integral = np.trapz(
        data_trapz_y,
        x = data_trapz_x
    )

    integral = -np.log(integral) * kBT
    return integral

def fes_1M_correct(
    volume: float, # angstrom**3
    temperature: float  # kelvin
) -> float:

    if volume == '1M':
        return 0

    Avogadro = 6.02214076e23
    volume_1M = 1e27/Avogadro
    kBT = T2kBT(temperature)
    #print(volume_1M/volume,'mol/L')
    correct = -kBT*math.log(volume_1M/volume)

    return correct

def get_pka(
    data_x,
    data_y,
    dict_coef: dict,
    tup_srange_tot: tuple,
    volume: float = None,
    temperature: float = None,
) -> (float, float):

    '''Compute pka.

    Args:
        dict_coef: 
            {
                (s0, s1): c0,
                (s1, s2): c1,
                ...
            }
        volume: volume
            V
        tup_srange_tot:
            (sl, sr)
        temperature: temperature
            T

    Returns:
          DeltaF 
        = 
          c0*[-kB*T*ln(P(s0,s1)/P(sl,sr))]
        + c1*[-kB*T*ln(P(s1,s2)/P(sl,sr))]
        + [-kB*T*ln(V_x/V)]

    '''

    h_tot = fes_integral(
        data_x, data_y,
        tup_xrange = tup_srange_tot,
        temperature = temperature,
    )
 
    deltag = 0
    for tup_srange, coef in dict_coef.items():
        h = fes_integral(
            data_x, data_y,
            tup_xrange = tup_srange,
            temperature = temperature,
        )
        deltag += coef*(h-h_tot)

    if not (volume is None):
        correct = fes_1M_correct(
            volume = volume,
            temperature = temperature
        )
        deltag += correct

    pka = deltag_to_pka(
        deltag = deltag,
        temperature = temperature,
    )

    return deltag, pka

def prob_to_deltag(
    prob: float,
    temperature: float, # kelvin
):

    kBT = T2kBT(temperature)
    return -kBT*np.log(prob)

def deltag_to_pka(
    deltag: float,
    temperature: float,
) -> float:

    kBT = T2kBT(temperature)
    return deltag/(kBT*math.log(10))

def pka_to_deltag(
    pka: float,
    temperature: float,
) -> float:

    kBT = T2kBT(temperature)
    return pka*(kBT*math.log(10))

def prob_to_pka(
    prob: float,
    temperature: float, # kelvin
    volume: float # angstrom**2
):
    deltag = prob_to_deltag(
        prob = prob,
        temperature = temperature
    )

    correct = fes_1M_correct(
        volume = volume,
        temperature = temperature
    )
    print('energy correct = ', correct)
    pka_correct = deltag_to_pka(
        deltag = correct,
        temperature = temperature,
    )
    print('pka correct = ', pka_correct)
    deltag += correct

    pka = deltag_to_pka(
        deltag = deltag,
        temperature = temperature,
    )

    return deltag, pka

def get_pka_time(
    dict_file: dict,
    dict_coef: dict,
    volume: float = None,
    tup_srange_tot: tuple = None,
    str_save: str = None,
    temperature: float = None,
):

    int_nfile = len(dict_file)
    np_pka = np.zeros(shape=(int_nfile), dtype=[('time', 'f4'), ('pka', 'f4')])
    np_deltag = np.zeros(shape=(int_nfile), dtype=[('time', 'f4'), ('deltag', 'f4')])
    print(temperature)
    for int_i,int_key in enumerate(dict_file):
        if not os.path.isfile(dict_file[int_key]):
            continue
        deltag, pka = get_pka(
            str_file = dict_file[int_key],
            dict_coef = dict_coef,
            tup_srange_tot = tup_srange_tot,
            temperature = temperature,
            volume = volume
        )
        np_pka[int_i] = (int_key, pka)
        np_deltag[int_i] = (int_key, deltag)

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
        print(str_file)
        with open(str_file, 'r') as file_open:
            str_line = file_open.readline()
            list_field = str_line.split()[2:]
        df_tmp = pd.read_csv(str_file, sep=' ', comment='#', names=list_field)
        df_data = pd.concat([df_data, df_tmp])

    return df_data

def colvar_plt(
    ax,
    header: str,
    list_data: list = None,
    timescale: float = 1.0,
    color: str = 'black',
    dict_color: dict = None,
    dict_label: dict = None,
    bool_scatter: bool = True,
    lw: float = None,
    scatters: float = None,
) -> None:

    if list_data is None:
        list_data = ['COLVAR']
    if dict_label is None:
        dict_label = {}

    df_data = colvar_reader(list_data)

    if bool_scatter:
        ax.scatter(df_data['time']*timescale, df_data[header], s=scatters, color=color, edgecolors='none', rasterized=True)
        if not(dict_color is None):
            for tup_range, color_i in dict_color.items():
                df_tmp = df_data.iloc[tup_range[0]:tup_range[1]]
                if not(tup_range in dict_label):
                    dict_label[tup_range] = None
                ax.scatter(df_tmp['time']*timescale, df_tmp[header], s=scatters, color=color_i, edgecolors='none', label=dict_label[tup_range], rasterized=True)
    else:
        ax.plot(df_data['time']*timescale, df_data[header], color=color, linewidth=lw)

    plot.set_lw( ax, lw )

def hills_sum(
    str_min: str,
    str_max: str,
    str_bin: str,
    str_cv: str = None,
    str_outfile: str = None,
    plm_in: str = None,
    plm_log: str = None,
    str_hills: str = 'HILLS'
):

    temperature, kBT = get_temperature(plm_in, plm_log)

    if not str_outfile:
        str_outfile = str_cv+'_fes.'

    str_cmd = f'source {os.environ["HOME"]}/.config/.tianff &&'
    str_cmd += 'source ${homedir}/.local/bin/bashrc_plm.sh ;'
    str_cmd += f'plumed sum_hills --hills {str_hills} --stride 20000 --outfile {str_outfile} --min {str_min} --max {str_max} --bin {str_bin} --negbias --kt {kBT}'
    if str_cv:
        str_cmd += f' --idw {str_cv}'

    subprocess_results = subprocess.run( str_cmd, shell=True, check=True, text=True, executable='/bin/bash')

    print(subprocess_results.stdout)

def T2kBT(
    temperature: float, # kelvin
) -> float:

    # Kcal/(mol*K)
    kB = 1.987204259e-3

    # Kcal/mol
    kBT = kB*temperature
    return kBT

def get_temperature(
    plm_in: str = '../plm.in',
    plm_log: str = '../plm.log'
) -> float:

    if plm_in:
        with open(plm_in, 'r') as fp:
            for str_line in fp:
                if 'TEMP=' in str_line:
                    list_line = str_line.split()
                    break
            for str_key in list_line:
                if 'TEMP=' == str_key[:5]:
                    temperature = int(str_key[5:])
                    break
        if not temperature:
            raise
        kBT = T2kBT(temperature)

    if plm_log:
        with open(plm_log, 'r') as fp:
            for str_line in fp:
                if 'kBT' in str_line:
                    list_line = str_line.split()
                    kBT_log = float(list_line[2])
                    break
        if not kBT_log:
            raise

    if plm_in and plm_log:
        if kBT-kBT_log > 1e4:
            raise

    if not (plm_in or plm_log):
        raise

    return temperature, kBT

