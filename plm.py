import math
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

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
    float_correct = -float_KbT*math.log(float_volume_1M/float_volume)

    return float_correct

def get_pka_x(
    str_file: str,
    tup_xrange1: tuple,
    tup_xrange2: tuple,
    float_volume: float,
    float_T: float = None,
    str_in: str = None,
    str_log: str = None
) -> (float, float):

    if not float_T:
        float_T, float_KbT = get_temperature( str_in, str_log )

    float_g1 = fes_integral(
        str_file = str_file,
        tup_xrange = tup_xrange1,
        float_T = float_T,
    )
    
    float_g2 = fes_integral(
        str_file = str_file,
        tup_xrange = tup_xrange2,
        float_T = float_T,
    )

    float_correct = fes_1M_correct(
        float_volume = float_volume,
        float_T = float_T
    )

    float_deltag = float_g2 - float_g1 + float_correct
    float_pka = deltag_to_pka(
        float_deltag = float_deltag,
        float_T = float_T,
    )

    return float_deltag, float_pka

def get_pka(
    str_file: str,
    dict_coef: dict,
    float_volume: float,
    tup_srange_tot: tuple,
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
    float_volume: float,
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
        str_save_pka = f'{str_save}_pka.csv'
        print(str_save_pka)
        np.savetxt(str_save_pka, np_pka, header=' '.join(np_pka.dtype.names))
        str_save_deltag = f'{str_save}_deltag.csv'
        print(str_save_deltag)
        np.savetxt(str_save_deltag, np_deltag, header=' '.join(np_deltag.dtype.names))

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

