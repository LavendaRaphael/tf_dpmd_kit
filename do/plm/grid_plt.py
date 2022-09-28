import numpy
from matplotlib import pyplot as plt
import numpy as np
import matplotlib

def grid_plt(
    list2d_file: list[list],
    str_xlabel: str,
    str_save: str=None,
    bool_minzero: bool=True,
    tup_ylim: tuple=None,
) -> None:

    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    fig, ax = plt.subplots()
    for list_file in list2d_file:
        str_file = list_file[0]
        str_label = list_file[1]

        np_data = np.loadtxt(str_file)
        if bool_minzero:
            ax.plot( np_data[:,0], np_data[:,1]-min(np_data[:,1]), label=str_label)
        else:
            ax.plot( np_data[:,0], np_data[:,1], label=str_label)
    
    ax.legend()
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel('FES (kJ/mol)')
    ax.set_xlim((1.1,7))
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')
#'''
grid_plt(
    list2d_file = [
        ['sum_hills/fes.dist_vp_c.4.dat', 'sum_hills'],
        ['reweight_bias/fes.dist_vp_c.grid', 'reweight_bias'],
        ['reweight_metad/fes.dist_vp_c.grid', 'reweight_metad'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.compare.pdf',
    tup_ylim = (None, 60)
)
#'''
'''
grid_plt(
    list2d_file = [
        ['fes.dist_vp_c.0.dat', '0.2ns'],
        ['fes.dist_vp_c.2.dat', '0.6ns'],
        ['fes.dist_vp_c.4.dat', '1.0ns'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.pdf',
)
#'''
'''
grid_plt(
    list2d_file = [
        ['analysis.0.fes.dist_vp_c.grid', '0.2ns'],
        ['analysis.2.fes.dist_vp_c.grid', '0.6ns'],
        ['fes.dist_vp_c.grid', '1.0ns'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.pdf',
)
#'''
plt.show()


