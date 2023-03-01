from matplotlib import pyplot as plt
import numpy as np
import matplotlib
from tf_dpmd_kit import plot

def time_ekin_plt(
    str_evp:  str='cp/cp.evp',
):

    with open(str_evp, 'r') as open_evp:
        list_header = open_evp.readline().split()[1:-2]
    np_evp = np.genfromtxt(str_evp, names=list_header)
    
    fig, ax = plt.subplots()
    ax.plot( np_evp['timeps'], np_evp['ekinc'], label = 'Electron')
    ax.plot( np_evp['timeps'], (np_evp['econs']-np_evp['etot']), label = 'Ion')
    
    ax.legend()
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel(f'Ekin ($E_h$)')
    ax.set_ylim(0, 1)

    return fig, ax

plot.set_rcparam()
fig, ax = time_ekin_plt()

cm = 1/2.54
plot.save(
    fig,
    str_save = 'time_ekin.svg',
    tup_size = (8.6*cm, 5*cm),
)

plt.show()

