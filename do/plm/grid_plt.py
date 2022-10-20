import numpy
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import rc
import os

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
        if not os.path.isfile(str_file):
            continue
        np_data = np.loadtxt(str_file)
        if bool_minzero:
            ax.plot( np_data[:,0], np_data[:,1]-min(np_data[:,1]), label=str_label, linewidth=2)
        else:
            ax.plot( np_data[:,0], np_data[:,1], label=str_label, linewidth=2)
    
    ax.legend()
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel('FES (kJ/mol)')
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.set_size_inches(6, 5)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)
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
    tup_xlim = (1.1, 8),
    tup_ylim = (-10, 70)
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
        ['sum_hills/fes.dist_vp_o_1_2.0.dat', '1ns'],
        ['sum_hills/fes.dist_vp_o_1_2.1.dat', '2ns'],
        ['sum_hills/fes.dist_vp_o_1_2.2.dat', '3ns'],
        ['sum_hills/fes.dist_vp_o_1_2.3.dat', '4ns'],
        ['sum_hills/fes.dist_vp_o_1_2.4.dat', '5ns'],
        ['sum_hills/fes.dist_vp_o_1_2.5.dat', '6ns'],
        ['sum_hills/fes.dist_vp_o_1_2.6.dat', '7ns'],
        ['sum_hills/fes.dist_vp_o_1_2.7.dat', '8ns'],
        ['sum_hills/fes.dist_vp_o_1_2.8.dat', '9ns'],
        ['sum_hills/fes.dist_vp_o_1_2.9.dat', '10ns'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'sum_hills/fes.dist_vp_c.nosft.pdf',
    bool_minzero = False,
    #tup_xlim = None
)
'''
#'''
grid_plt(
    list2d_file = [
        ['sum_hills/fes.dist_vp_o_1_2.0.dat', '1ns'],
        ['sum_hills/fes.dist_vp_o_1_2.1.dat', '2ns'],
        ['sum_hills/fes.dist_vp_o_1_2.2.dat', '3ns'],
        ['sum_hills/fes.dist_vp_o_1_2.3.dat', '4ns'],
        ['sum_hills/fes.dist_vp_o_1_2.4.dat', '5ns'],
        ['sum_hills/fes.dist_vp_o_1_2.5.dat', '6ns'],
        ['sum_hills/fes.dist_vp_o_1_2.6.dat', '7ns'],
        ['sum_hills/fes.dist_vp_o_1_2.7.dat', '8ns'],
        ['sum_hills/fes.dist_vp_o_1_2.8.dat', '9ns'],
        ['sum_hills/fes.dist_vp_o_1_2.9.dat', '10ns'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_save = 'sum_hills/fes.dist_vp_c.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 70)
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
    tup_xlim = (1.1, 8),
    tup_ylim = (None, 70)
)
#'''
plt.show()


