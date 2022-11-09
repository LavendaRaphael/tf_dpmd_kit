import numpy
from matplotlib import pyplot as plt
import numpy as np
import matplotlib
import pandas as pd
import seaborn as sns
import matplotlib.ticker as plticker


def grid2d_plt(
    str_file: str,
    str_xlabel: str = None,
    str_ylabel: str = None,
    bool_minus: bool = False,
    bool_minzero: bool = False,
    str_save: str = None
) -> None:

    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    with open(str_file, 'r') as file_obj:
        list_field = file_obj.readline().split()[2:]
        while True:
            list_line = file_obj.readline().split()
            if list_line[0] != '#!':
                break
            if list_line[2] == f'nbins_{list_field[0]}':
                int_nbin_x = int(list_line[3])
            elif list_line[2] == f'nbins_{list_field[1]}':
                int_nbin_y = int(list_line[3])
        if not (int_nbin_x and int_nbin_y):
            raise

    np_data = np.loadtxt(str_file)
    x = np_data[:,0]
    y = np_data[:,1]
    z = np_data[:,2]
    if bool_minus:
        z = -z
    if bool_minzero:
        z -= min(z)
    z = z.reshape(int_nbin_y, int_nbin_x)

    fig, ax = plt.subplots()
    image = ax.imshow(
        z,
        origin = 'lower',
        extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
        #cmap='coolwarm',
        cmap='Reds',
        aspect='auto'
    )
    fig.colorbar(image)
    #loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
    #ax.xaxis.set_major_locator(loc)

    if not str_xlabel:
        str_xlabel = list_field[0]
    if not str_ylabel:
        str_ylabel = list_field[1]
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel(str_ylabel)
    #fig.set_size_inches(8, 8)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

grid2d_plt(
    str_file = 'cv_2d.5.dat',
    bool_minus = True,
    bool_minzero = True,
    #str_xlabel = r'dh(O$_0$H$_0$) (rad)',
    #str_ylabel = r'dh(O$_1$V$_H$) (rad)',
    str_save = 'cv_2d.5.pdf'
)

plt.show()


