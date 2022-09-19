import numpy
from matplotlib import pyplot as plt
import numpy as np
import matplotlib
import pandas as pd
import seaborn as sns
import matplotlib.ticker as plticker


def grid2d_plt(
    str_file: str,
    str_xlabel: str,
    str_ylabel: str,
    str_save: str=None
) -> None:

    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    with open(str_file, 'r') as file_obj:
        list_field = file_obj.readline().split()[2:]
        for i in range(3):
            file_obj.readline()
        list_line = file_obj.readline().split()[2:]
        if list_line[0] == f'nbins_{list_field[0]}':
            int_nbin_x = int(list_line[1]) + 1
        else:
            raise
        for i in range(3):
            file_obj.readline()
        list_line = file_obj.readline().split()[2:]
        if list_line[0] == f'nbins_{list_field[1]}':
            int_nbin_y = int(list_line[1]) + 1
        else:
            raise

    np_data = np.loadtxt(str_file)
    x = np_data[:,0]
    y = np_data[:,1]
    z = np_data[:,2]
    z = z.reshape(int_nbin_y, int_nbin_x)

    fig, ax = plt.subplots()
    image = ax.imshow(
        z,
        origin = 'lower',
        extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
        cmap='coolwarm'
    )
    fig.colorbar(image)
    loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)

    ax.set_xlabel(str_xlabel)
    ax.set_ylabel(str_ylabel)
    #fig.set_size_inches(8, 4)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

grid2d_plt(
    str_file = 'fes.dh_o_h.grid',
    str_xlabel = r'dh(O$_0$H$_0$) (rad)',
    str_ylabel = r'dh(O$_1$V$_H$) (rad)',
    str_save = 'fes.dh_o_h.pdf'
)

plt.show()


