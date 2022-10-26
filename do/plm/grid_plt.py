import numpy
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import rc
import os

def grid_plt(
    list2d_file: list[list],
    str_xlabel: str,
    str_ylabel: str = 'FES (kJ/mol)',
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
            ax.plot( np_data[:,0], np_data[:,1]-min(np_data[:,1]), label=str_label, linewidth=1)
        else:
            ax.plot( np_data[:,0], np_data[:,1], label=str_label, linewidth=1)
    
    ax.legend()
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel(str_ylabel)
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.set_size_inches(6, 5)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)
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
'''
grid_plt(
    list2d_file = [
        #['fes.dist_vp_o_1_2.0.dat', '1ns'],
        #['fes.dist_vp_o_1_2.1.dat', '2ns'],
        #['fes.dist_vp_o_1_2.2.dat', '3ns'],
        #['fes.dist_vp_o_1_2.3.dat', '4ns'],
        #['fes.dist_vp_o_1_2.4.dat', '5ns'],
        ['fes.dist_vp_o_1_2.5.dat', '6ns'],
        #['fes.dist_vp_o_1_2.6.dat', '7ns'],
        ['fes.dist_vp_o_1_2.7.dat', '8ns'],
        #['fes.dist_vp_o_1_2.8.dat', '9ns'],
        ['fes.dist_vp_o_1_2.9.dat', '10ns'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = 'fes.dist_vp_c.nosft.pdf',
    bool_minzero = False,
    #tup_xlim = None
)
#'''
'''
grid_plt(
    list2d_file = [
        #['fes.dist_vp_o_1_2.0.dat', '1ns'],
        #['fes.dist_vp_o_1_2.1.dat', '2ns'],
        #['fes.dist_vp_o_1_2.2.dat', '3ns'],
        #['fes.dist_vp_o_1_2.3.dat', '4ns'],
        #['fes.dist_vp_o_1_2.4.dat', '5ns'],
        ['fes.dist_vp_o_1_2.5.dat', '6ns'],
        #['fes.dist_vp_o_1_2.6.dat', '7ns'],
        ['fes.dist_vp_o_1_2.7.dat', '8ns'],
        #['fes.dist_vp_o_1_2.8.dat', '9ns'],
        ['fes.dist_vp_o_1_2.9.dat', '10ns'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = 'fes.dist_vp_c.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 70)
)
#'''
'''
str_tmp = 'dist_vp_c_hgx'
grid_plt(
    list2d_file = [
        #[f'analysis.0.{str_tmp}.grid', '1ns'],
        #[f'analysis.1.{str_tmp}.grid', '2ns'],
        #[f'analysis.2.{str_tmp}.grid', '3ns'],
        #[f'analysis.3.{str_tmp}.grid', '4ns'],
        #[f'analysis.4.{str_tmp}.grid', '5ns'],
        [f'analysis.5.{str_tmp}.grid', '6ns'],
        #[f'analysis.6.{str_tmp}.grid', '7ns'],
        [f'analysis.7.{str_tmp}.grid', '8ns'],
        #[f'analysis.8.{str_tmp}.grid', '9ns'],
        [f'{str_tmp}.grid', '10ns'],
    ],
    str_xlabel = r'R(V$_P$C) (Å)',
    str_ylabel = 'Probility Density',
    str_save = f'{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 70)
)
#'''
'''
str_tmp = 'dist_vp_c_fes'
grid_plt(
    list2d_file = [
        #[f'analysis.0.{str_tmp}.grid', '1ns'],
        #[f'analysis.1.{str_tmp}.grid', '2ns'],
        #[f'analysis.2.{str_tmp}.grid', '3ns'],
        #[f'analysis.3.{str_tmp}.grid', '4ns'],
        #[f'analysis.4.{str_tmp}.grid', '5ns'],
        [f'analysis.5.{str_tmp}.grid', '6ns'],
        #[f'analysis.6.{str_tmp}.grid', '7ns'],
        [f'analysis.7.{str_tmp}.grid', '8ns'],
        #[f'analysis.8.{str_tmp}.grid', '9ns'],
        [f'{str_tmp}.grid', '10ns'],
    ],
    str_xlabel = r'R(V$_P$C) (Å)',
    str_save = f'{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 70)
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes'
grid_plt(
    list2d_file = [
        #[f'analysis.0.{str_tmp}.grid', '1ns'],
        #[f'analysis.1.{str_tmp}.grid', '2ns'],
        #[f'analysis.2.{str_tmp}.grid', '3ns'],
        #[f'analysis.3.{str_tmp}.grid', '4ns'],
        #[f'analysis.4.{str_tmp}.grid', '5ns'],
        [f'analysis.5.{str_tmp}.grid', '6ns'],
        #[f'analysis.6.{str_tmp}.grid', '7ns'],
        [f'analysis.7.{str_tmp}.grid', '8ns'],
        #[f'analysis.8.{str_tmp}.grid', '9ns'],
        [f'{str_tmp}.grid', '10ns'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = f'{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 70)
)
#'''
plt.show()


