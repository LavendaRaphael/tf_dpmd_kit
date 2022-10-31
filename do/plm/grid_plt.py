import numpy
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import rc
import os
from matplotlib import cm

def grid_plt(
    list2d_file: list[list],
    str_xlabel: str,
    str_ylabel: str,
    str_save: str=None,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    tup_colormap: tuple = None,
    bool_minzero: bool = False,
    bool_maxzero: bool = False,
    bool_minus: bool = False,
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    if tup_colormap:
        sm = plt.cm.ScalarMappable(cmap=cm.coolwarm, norm=plt.Normalize(vmin=tup_colormap[0], vmax=tup_colormap[1]))

    fig, ax = plt.subplots()
    for list_file in list2d_file:
        str_file = list_file[0]
        str_label = list_file[1]
        if not os.path.isfile(str_file):
            continue

        if tup_colormap:
            color = sm.to_rgba(int(str_label[:-1]))
        else:
            color = None

        np_data = np.loadtxt(str_file)
        np_data_y = np_data[:,1]
        if bool_minus:
            np_data_y *= -1

        if bool_minzero:
            np_data_y -= min(np_data_y)

        if bool_maxzero:
            np_data_y -= max(np_data_y)

        ax.plot( np_data[:,0], np_data_y, label=str_label, linewidth=2, color=color)
    
    ax.legend()
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel(str_ylabel)
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.set_size_inches(9, 7)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)
#------------------------------------------------------------------[Temperature]
#'''
str_tmp = 'dist_vp_c_deltag'
grid_plt(
    list2d_file = [
        [f'290K/reweight_bias/{str_tmp}.csv', '290K'],
        [f'310K/reweight_bias/{str_tmp}.csv', '310K'],
        [f'330K/reweight_bias/{str_tmp}.csv', '330K'],
        [f'350K/reweight_bias/{str_tmp}.csv', '350K'],
        [f'370K/reweight_bias/{str_tmp}.csv', '370K'],
        [f'390K/reweight_bias/{str_tmp}.csv', '390K'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = r'$\Delta$G (kJ/mol)',
    str_save = f'reweight_bias.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (-10, 70),
    tup_colormap = (290, 390)
)
#'''
'''
str_tmp = 'dist_vp_c_pka'
grid_plt(
    list2d_file = [
        [f'290K/reweight_bias/{str_tmp}.csv', '290K'],
        [f'310K/reweight_bias/{str_tmp}.csv', '310K'],
        [f'330K/reweight_bias/{str_tmp}.csv', '330K'],
        [f'350K/reweight_bias/{str_tmp}.csv', '350K'],
        [f'370K/reweight_bias/{str_tmp}.csv', '370K'],
        [f'390K/reweight_bias/{str_tmp}.csv', '390K'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = 'pKa',
    str_save = f'reweight_bias.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    tup_ylim = (None, 4.2),
    tup_colormap = (290, 390)
)
#'''
'''
str_tmp='fes.dist_vp_c'
grid_plt(
    list2d_file = [
        [f'290K/sum_hills/{str_tmp}.grid', '290K'],
        [f'310K/sum_hills/{str_tmp}.grid', '310K'],
        [f'330K/sum_hills/{str_tmp}.grid', '330K'],
        [f'350K/sum_hills/{str_tmp}.grid', '350K'],
        [f'370K/sum_hills/{str_tmp}.grid', '370K'],
        [f'390K/sum_hills/{str_tmp}.grid', '390K'],
    ],
    str_xlabel = r'R(CV$_P$) (Å)',
    str_ylabel = 'Hills (kJ/mol)',
    str_save = f'sum_hills.{str_tmp}.pdf',
    bool_minzero = True,
    bool_minus = True,
    tup_xlim = (0,10),
    tup_ylim = (None, None),
    tup_colormap = (290, 390)
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes'
grid_plt(
    list2d_file = [
        [f'290K/reweight_bias/{str_tmp}.grid', '290K'],
        [f'310K/reweight_bias/{str_tmp}.grid', '310K'],
        [f'330K/reweight_bias/{str_tmp}.grid', '330K'],
        [f'350K/reweight_bias/{str_tmp}.grid', '350K'],
        [f'370K/reweight_bias/{str_tmp}.grid', '370K'],
        [f'390K/reweight_bias/{str_tmp}.grid', '390K'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_ylabel = 'FES (kJ/mol)',
    str_save = f'reweight_bias.{str_tmp}.pdf',
    bool_minzero = True,
    tup_xlim = (0, 14),
    tup_ylim = (None, 60),
    tup_colormap = (290, 390),
)
#'''
'''
str_tmp = 'dist_vp_c_fes'
#str_tmp = 'fes.dist_vp_c'
grid_plt(
    list2d_file = [
        [f'290K/reweight_bias/{str_tmp}.grid', '290K'],
        [f'310K/reweight_bias/{str_tmp}.grid', '310K'],
        [f'330K/reweight_bias/{str_tmp}.grid', '330K'],
        [f'350K/reweight_bias/{str_tmp}.grid', '350K'],
        [f'370K/reweight_bias/{str_tmp}.grid', '370K'],
        [f'390K/reweight_bias/{str_tmp}.grid', '390K'],
    ],
    str_xlabel = r'R(V$_P$C) (Å)',
    str_ylabel = 'FES (kJ/mol)',
    str_save = f'reweight_bias.{str_tmp}.pdf',
    tup_xlim = (0, 10),
    tup_ylim = (None, 60),
    bool_minzero = True,
    tup_colormap = (290, 390)
)
#'''
#----------------------------------------------------------[Compare]

'''
str_tmp = 'dist_vp_c_pka'
grid_plt(
    list2d_file = [
        [f'../../../03.390K_hills_reweight/390K/reweight_bias/{str_tmp}.csv', 'Biasfactor 5.0'],
        [f'{str_tmp}.csv', 'Biasfactor 10.0'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = 'pKa',
    str_save = f'{str_tmp}.compare.pdf'
)
#'''
'''
str_tmp = 'dist_vp_c_deltag'
grid_plt(
    list2d_file = [
        [f'../../../03.390K_hills_reweight/390K/reweight_bias/{str_tmp}.csv', 'Biasfactor 5.0'],
        [f'{str_tmp}.csv', 'Biasfactor 10.0'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = r'$\Delta$G (kJ/mol)',
    str_save = f'{str_tmp}.compare.pdf'
)
#'''
'''
str_tmp = 'dist_vp_c_fes'
grid_plt(
    list2d_file = [
        [f'../../../03.390K_hills_reweight/390K/reweight_bias/{str_tmp}.grid', 'Biasfactor 5.0'],
        [f'{str_tmp}.grid', 'Biasfactor 10.0'],
    ],
    str_xlabel = r'R(V$_P$C) (Å)',
    str_save = f'{str_tmp}.compare.pdf',
    bool_minzero = True,
    str_ylabel = 'FES (kJ/mol)',
    tup_xlim = (0,10),
    tup_ylim = (None, 80),
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes'
grid_plt(
    list2d_file = [
        [f'../../../03.390K_hills_reweight/390K/reweight_bias/{str_tmp}.grid', 'Biasfactor 5.0'],
        [f'{str_tmp}.grid', 'Biasfactor 10.0'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = f'{str_tmp}.compare.pdf',
    bool_minzero = True,
    str_ylabel = 'FES (kJ/mol)',
    tup_xlim = (0,14),
    tup_ylim = (None, 80),
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes'
grid_plt(
    list2d_file = [
        ['../../../02.390K_hills/390K/sum_hills/fes.dist_vp_o_1_2.10.dat', 'Biasfactor 5.0'],
        [f'{str_tmp}.5.dat', 'Biasfactor 10.0'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = f'{str_tmp}.compare.pdf',
    bool_minus = True,
    bool_maxzero = True,
    str_ylabel = 'Hills (kJ/mol)',
    tup_xlim = (0,14)
)
#'''
#----------------------------------------------------------[Sum hills]

'''
str_tmp = 'dist_vp_o_1_2_fes'
#str_tmp = 'fes.dist_vp_o_1_2'
grid_plt(
    list2d_file = [
        [f'{str_tmp}.0.dat', '1ns'],
        [f'{str_tmp}.1.dat', '2ns'],
        [f'{str_tmp}.2.dat', '3ns'],
        [f'{str_tmp}.3.dat', '4ns'],
        [f'{str_tmp}.4.dat', '5ns'],
        [f'{str_tmp}.5.dat', '6ns'],
        #[f'{str_tmp}.6.dat', '7ns'],
        [f'{str_tmp}.7.dat', '8ns'],
        #[f'{str_tmp}.8.dat', '9ns'],
        [f'{str_tmp}.9.dat', '10ns'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = f'{str_tmp}.nosft.pdf',
    bool_minus = True,
    bool_minzero = True,
    str_ylabel = 'Hills (kJ/mol)',
    tup_xlim = (0,14)
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes'
#str_tmp = 'fes.dist_vp_o_1_2'
grid_plt(
    list2d_file = [
        [f'{str_tmp}.0.dat', '1ns'],
        [f'{str_tmp}.1.dat', '2ns'],
        [f'{str_tmp}.2.dat', '3ns'],
        [f'{str_tmp}.3.dat', '4ns'],
        [f'{str_tmp}.4.dat', '5ns'],
        [f'{str_tmp}.5.dat', '6ns'],
        #[f'{str_tmp}.6.dat', '7ns'],
        [f'{str_tmp}.7.dat', '8ns'],
        #[f'{str_tmp}.8.dat', '9ns'],
        [f'{str_tmp}.9.dat', '10ns'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = f'{str_tmp}.pdf',
    bool_minus = True,
    bool_maxzero = True,
    str_ylabel = 'Hills (kJ/mol)',
    tup_xlim = (0, 14),
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
