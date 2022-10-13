import numpy
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import rc

def grid_plt(
    list2d_file: list[list],
    str_xlabel: str,
    str_save: str=None,
    bool_minzero: bool=True,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

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
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')
'''
grid_plt(
    list2d_file = [
        ['290K/fes/sum_hills/fes.dist_vp_c.grid', '290K'],
        ['310K/fes/sum_hills/fes.dist_vp_c.grid', '310K'],
        ['330K/fes/sum_hills/fes.dist_vp_c.grid', '330K'],
        ['350K/fes/sum_hills/fes.dist_vp_c.grid', '350K'],
        ['370K/fes/sum_hills/fes.dist_vp_c.grid', '370K'],
        ['390K/fes/sum_hills/fes.dist_vp_c.grid', '390K'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.sum_hills.pdf',
    tup_ylim = (None, 70)
)
#'''
'''
grid_plt(
    list2d_file = [
        ['290K/fes/reweight_bias/fes.dist_vp_c.grid', '290K'],
        ['310K/fes/reweight_bias/fes.dist_vp_c.grid', '310K'],
        ['330K/fes/reweight_bias/fes.dist_vp_c.grid', '330K'],
        ['350K/fes/reweight_bias/fes.dist_vp_c.grid', '350K'],
        ['370K/fes/reweight_bias/fes.dist_vp_c.grid', '370K'],
        ['390K/fes/reweight_bias/fes.dist_vp_c.grid', '390K'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.reweight_bias.pdf',
    tup_ylim = (None, 70)
)
#'''
'''
grid_plt(
    list2d_file = [
        ['sum_hills/fes.dist_vp_c.10.dat', 'sum_hills'],
        ['reweight_bias/fes.dist_vp_c.grid', 'reweight_bias'],
        #['reweight_metad/fes.dist_vp_c.grid', 'reweight_metad'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'fes.dist_vp_c.compare.pdf',
    tup_ylim = (None, 70)
)
#'''
#'''
grid_plt(
    list2d_file = [
        ['sum_hills/fes.dist_vp_c.1.dat', '2ns'],
        ['sum_hills/fes.dist_vp_c.3.dat', '4ns'],
        ['sum_hills/fes.dist_vp_c.5.dat', '6ns'],
        ['sum_hills/fes.dist_vp_c.7.dat', '8ns'],
        ['sum_hills/fes.dist_vp_c.9.dat', '10ns']
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'sum_hills/fes.dist_vp_c.minzero.pdf',
    bool_minzero = False,
    tup_xlim = (0, 10)
)
#'''
'''
grid_plt(
    list2d_file = [
        ['sum_hills/fes.dist_vp_c.1.dat', '2ns'],
        ['sum_hills/fes.dist_vp_c.3.dat', '4ns'],
        ['sum_hills/fes.dist_vp_c.5.dat', '6ns'],
        ['sum_hills/fes.dist_vp_c.7.dat', '8ns'],
        ['sum_hills/fes.dist_vp_c.9.dat', '10ns']
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'sum_hills/fes.dist_vp_c.pdf',
    tup_xlim = (0.85, 8)
    tup_ylim = (None, 70)
)
#'''
'''
grid_plt(
    list2d_file = [
        ['reweight_bias/analysis.1.fes.dist_vp_c.grid', '2ns'],
        ['reweight_bias/analysis.3.fes.dist_vp_c.grid', '4ns'],
        ['reweight_bias/analysis.5.fes.dist_vp_c.grid', '6ns'],
        ['reweight_bias/analysis.7.fes.dist_vp_c.grid', '8ns'],
        ['reweight_bias/fes.dist_vp_c.grid', '10ns'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'reweight_bias/fes.dist_vp_c.pdf',
    tup_ylim = (None, 70)
)
#'''
plt.show()


