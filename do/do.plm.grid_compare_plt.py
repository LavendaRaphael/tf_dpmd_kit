import numpy
from matplotlib import pyplot as plt
import numpy as np
import matplotlib


def grid_plt(
    list2d_file: list[list],
    str_xlabel: str,
    str_save: str=None
) -> None:

    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    fig, ax = plt.subplots()
    for list_file in list2d_file:
        str_file = list_file[0]
        str_label = list_file[1]

        np_data = np.loadtxt(str_file)
        ax.plot( np_data[:,0], np_data[:,1]-(np_data[34,1]), label=str_label)
    
    ax.legend()
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel('FES (kJ/mol)')
    ax.set_xlim((1.1,7))
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

grid_plt(
    list2d_file = [
        ['sum_hills/fes.dist_vp_c.4.dat', 'sum_hills'],
        ['reweight_bias/fes.dist_vp_c.grid', 'reweight_bias'],
        ['reweight_metad/fes.dist_vp_c.grid', 'reweight_metad'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.compare.pdf',
)

plt.show()


