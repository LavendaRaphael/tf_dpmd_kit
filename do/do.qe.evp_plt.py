from matplotlib import pyplot
import numpy as np
import matplotlib

def def_time_ekin_plt(
    str_evp:  str='cp/cp.evp',
    str_save: str='time_ekin.pdf'
) -> None:

    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]
    
    with open(str_evp, 'r') as open_evp:
        list_header = open_evp.readline().split()[1:-2]
    np_evp = np.genfromtxt(str_evp, names=list_header)
    
    fig, ax = pyplot.subplots()
    ax.plot( np_evp['timeps'], np_evp['ekinc'], label = 'Ekin_e', linewidth=2)
    ax.plot( np_evp['timeps'], (np_evp['econs']-np_evp['etot'])/10, label = 'Ekin_ion/10', linewidth=2)
    
    ax.legend()
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel(f'Ekin ($E_h$)')
    ax.set_ylim(0, 0.1)
    fig.set_size_inches(8, 4)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

def_time_ekin_plt()

pyplot.show()

