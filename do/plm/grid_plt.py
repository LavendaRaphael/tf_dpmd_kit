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
    ax.set_xlim((0.85,8))
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')
'''
grid_plt(
    list2d_file = [
        ['sum_hills/fes.dist_vp_c.10.dat', 'sum_hills'],
        ['reweight_bias/fes.dist_vp_c.grid', 'reweight_bias'],
        ['reweight_metad/fes.dist_vp_c.grid', 'reweight_metad'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.compare.pdf',
    tup_ylim = (None, 70)
)
#'''
'''
grid_plt(
    list2d_file = [
        ['fes.dist_vp_c.1.dat', '2ns'],
        ['fes.dist_vp_c.3.dat', '4ns'],
        ['fes.dist_vp_c.5.dat', '6ns'],
        ['fes.dist_vp_c.7.dat', '8ns'],
        ['fes.dist_vp_c.9.dat', '10ns']
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.pdf',
    tup_ylim = (None, 70)
)
#'''
#'''
grid_plt(
    list2d_file = [
        #['analysis.0.fes.dist_vp_c.grid', '2ns'],
        ['analysis.1.fes.dist_vp_c.grid', '2ns'],
        #['analysis.2.fes.dist_vp_c.grid', '2ns'],
        ['analysis.3.fes.dist_vp_c.grid', '4ns'],
        #['analysis.4.fes.dist_vp_c.grid', '4ns'],
        ['analysis.5.fes.dist_vp_c.grid', '6ns'],
        #['analysis.6.fes.dist_vp_c.grid', '6ns'],
        ['analysis.7.fes.dist_vp_c.grid', '8ns'],
        #['analysis.8.fes.dist_vp_c.grid', '8ns'],
        ['fes.dist_vp_c.grid', '10ns'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.pdf',
    tup_ylim = (None, 70)
)
#'''
plt.show()


