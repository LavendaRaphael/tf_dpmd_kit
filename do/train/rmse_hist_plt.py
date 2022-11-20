import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import json

def def_plt(
    str_file: str,
    int_natoms: int = None,
    str_save: str = None,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    int_bins: int = 'auto',
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})
 
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

    float_rmse_e, float_rmse_f, float_rmse_v = get_rmse()

    if (str_mode=='e'):
        # per atom
        np_data /= int_natoms
        # eV to meV
        np_data *= 1000
        float_rmse_e *= 1000
        del_data = np_data[:,1] - np_data[:,0]
        ax.set_xlabel( r'E$_{DP}$-E$_{DFT}$ (meV/atom)')
        str_label = f'Energy RMSE = {float_rmse_e:.3f} meV/atom'

    elif (str_mode=='f'):
        del_data_xyz = np_data[:,3:6] - np_data[:,0:3]
        del_data = np.linalg.norm( del_data_xyz, axis=1 )
        ax.set_xlabel(r'|F$_{DP}$-F$_{DFT}$| (eV/Å)')
        str_label = f'Force RMSE = {float_rmse_f:.3f} eV/Å'

    ax.hist(
        del_data,
        label = str_label,
        bins = int_bins,
        density = True
    )

    ax.legend()
    ax.set_ylabel('Probabilty Density')
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

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

#'''
def_plt(
    str_file = 'dptest.e.out',
    int_natoms = 384,
    tup_xlim = (-1.2,1.2),
    tup_ylim = (0,1.6),
    str_save = 'dptest.e.pdf',
)
#'''
'''
def_plt(
    str_file = 'dptest.f.out',
    tup_xlim = (0,0.4),
    #tup_ylim = (0,10),
    str_save = 'dptest.f.pdf',
)
#'''

plt.show()
