import numpy as np
from matplotlib import pyplot as plt
import matplotlib

def def_plt(
    float_rmse: float,
    str_file: str,
    int_natoms: int=None,
    str_save: str=None,
    tup_xlim: tuple=None,
    int_bins: int='auto',
) -> None:
    
    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    fig, ax = plt.subplots()

    with open(str_file, 'r') as file_open:
        list_line = file_open.readline().split()
        if list_line[-1] == 'pred_e':
            str_mode = 'e'
            if not int_natoms:
                raise
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
        str_xlabel = r'E$_{DP}$-E$_{DFT}$ (meV/atom)'
        str_label = f'Energy RMSE = {float_rmse:.3f} meV/atom'
    elif (str_mode=='f'):
        del_data_xyz = np_data[:,3:6] - np_data[:,0:3]
        del_data = np.linalg.norm( del_data_xyz, axis=1 )
        del_data *= 1000
        str_xlabel = r'|F$_{DP}$-F$_{DFT}$| (meV/Å)'
        str_label = f'Force RMSE = {float_rmse:.1f} meV/Å'

    float_std = np.std(del_data)
    float_mean = np.mean(del_data) 
    ax.hist(
        del_data,
        label = str_label,
        bins = int_bins,
        density = True
    )

    ax.legend()
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel('Probabilty Density')
    ax.set_xlim(tup_xlim)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

#'''
def_plt(
    str_file = 'dptest.e.out',
    int_natoms = 384,
    float_rmse = 0.440,
    tup_xlim = (-1.5,1.5),
    str_save = 'dptest.e.pdf'
)
#'''
#'''
def_plt(
    str_file = 'dptest.f.out',
    float_rmse = 60.0,
    tup_xlim = (0,250),
    str_save = 'dptest.f.pdf'
)
#'''

plt.show()
