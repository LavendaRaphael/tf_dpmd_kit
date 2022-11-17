import numpy as np
from matplotlib import pyplot as plt
import matplotlib

def def_plt(
    float_rmse: float,
    str_file: str,
    int_natoms: int=None,
    str_save: str=None,
) -> None:
    
    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

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
        np_data -= np.average( np_data[:,0] )

        fig, ax = plt.subplots()
        ax.set_xlabel('DFT energy (meV/atom)')
        ax.set_ylabel('DNN energy (meV/atom)')
        ax.plot(
            np_data[:,0],
            np_data[:,1],
            label = f'Energy RMSE = {float_rmse:.3f} meV/atom',
            marker='o',
            linestyle = '',
            markersize=2,
        )
    elif (str_mode=='f'):
        np_data *= 1000

        fig, ax = plt.subplots()
        for in []:
            ax.plot(
                np_data[:,0],
                np_data[:,1],
                label = f'Energy RMSE = {float_rmse:.3f} meV/atom',
                marker='o',
                linestyle = '',
                markersize=2,
            )
        str_xlabel = 'DFT force (meV/Å)'
        str_ylabel = 'DNN force (meV/Å)'
        str_label = f'Force RMSE = {float_rmse:.1f} meV/Å'

    ax.axline([0, 0], [1, 1], color='black', linestyle='--')
    ax.legend()
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

#'''
def_plt(
    str_file = 'dptest.e.out',
    int_natoms = 384,
    float_rmse = 0.440,
)
#'''
'''
def_plt(
    str_file = 'dptest.f.out',
    float_rmse = 60.0,
    str_save = 'dptest.f.pdf'
)
#'''

plt.show()
