import numpy as np
from matplotlib import pyplot as plt
import matplotlib

def def_plt(
    list_file: list,
    str_save: str=None
) -> None:
    
    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    fig, ax = plt.subplots()

    for tup_file in list_file:
        np_data = np.loadtxt(tup_file[0])
        np_id = np.arange(len(np_data))
        del_data_xyz = np_data[:,3:6] - np_data[:,0:3]
        del_data = np.linalg.norm( del_data_xyz, axis=1 )
        del_data *= 1000

        ax.scatter(np_id, del_data, label=tup_file[1])

    ax.legend()
    ax.set_ylabel(r'|F$_{DP}$-F$_{DFT}$| (meV/Å)')
    ax.set_xlabel('Atom index')

    # filter large
    #for int_i in range(len(del_data)):
    #    if del_data[int_i] > 300:
    #        print(int_i, del_data[int_i])

    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

#'''
def_plt(
    list_file = [
        ('dptest.f.out', ''),
    ],
)
#'''
'''
def_plt(
    list_file = [
        ('snap_361413/dptest/dptest.f.out', 'snap_361413'),
        ('snap_362000/dptest/dptest.f.out', 'snap_362000'),
        ('snap_370000/dptest/dptest.f.out', 'snap_370000'),
        ('snap_399700/dptest/dptest.f.out', 'snap_399700'),
    ],
    str_save = 'dptest.f.pdf'
)
#'''

plt.show()
