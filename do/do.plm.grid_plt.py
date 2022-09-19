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
        ax.plot( np_data[:,0], np_data[:,1], label=str_label)
    
    ax.legend()
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel('FES (kJ/mol)')
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

grid_plt(
    list2d_file = [
        ['fes.dist_vp_c.0.dat', '0.2ns'],
        ['fes.dist_vp_c.1.dat', '0.4ns'],
        ['fes.dist_vp_c.2.dat', '0.6ns'],
        ['fes.dist_vp_c.3.dat', '0.8ns'],
        ['fes.dist_vp_c.4.dat', '1.0ns'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.pdf',
)
'''
grid_plt(
    list2d_file = [
        ['analysis.0.fes.dist_vp_c.grid', '0.2ns'],
        ['analysis.1.fes.dist_vp_c.grid', '0.4ns'],
        ['analysis.2.fes.dist_vp_c.grid', '0.6ns'],
        ['analysis.3.fes.dist_vp_c.grid', '0.8ns'],
        ['fes.dist_vp_c.grid', '1.0ns'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.pdf',
)
#'''

plt.show()


