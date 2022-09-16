import numpy
from matplotlib import pyplot as plt
import numpy as np
import matplotlib
import pandas as pd
import seaborn as sns
import matplotlib.ticker as plticker

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

def grid2d_plt(
    str_file: str,
    str_save: str=None
):
    fig, ax = plt.subplots()

    with open(str_file, 'r') as open_file:
        list_header = open_file.readline().split()[2:]
    np_data = np.genfromtxt(str_file, names=list_header)
    
    x = np_data['dhx_o_0_h']
    y = np_data['dhx_o_1_2_vh']
    z = np_data['fes_dh_o_h']

    N = int(len(z)**.5)
    z = np.flip(z.reshape(N, N), axis=0)
    image = ax.imshow(z, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
            cmap='coolwarm')
    fig.colorbar(image)
    loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)

    ax.set_xlabel(r'dh (O$_0$H$_0$)')
    ax.set_ylabel(r'dh (O$_1$V$_H$)')
    #fig.set_size_inches(8, 4)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

grid2d_plt(
    str_file = 'fes.dh_o_h.grid',
    str_save = 'fes.dh_o_h.pdf'
)

plt.show()


