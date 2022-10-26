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
    bool_minizero: bool = False
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
        np_data = np.loadtxt(str_file)
        if tup_colormap:
            color = sm.to_rgba(int(str_label[:-1]))
        else:
            color = None
        if bool_minizero:
            ax.plot( np_data[:,0], np_data[:,1]-min(np_data[:,1]), label=str_label, linewidth=1, color=color)
        else:
            ax.plot( np_data[:,0], np_data[:,1], label=str_label, linewidth=2, color=color)
    
    ax.legend()
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel(str_ylabel)
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.set_size_inches(9, 7)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)
#'''
str_tmp = 'dist_vp_c_deltag'
grid_plt(
    list2d_file = [
        [f'290K.0/reweight_bias/{str_tmp}.csv', '290K'],
        [f'310K.0/reweight_bias/{str_tmp}.csv', '310K'],
        [f'330K.0/reweight_bias/{str_tmp}.csv', '330K'],
        [f'350K.0/reweight_bias/{str_tmp}.csv', '350K'],
        [f'370K.0/reweight_bias/{str_tmp}.csv', '370K'],
        [f'390K.0/reweight_bias/{str_tmp}.csv', '390K'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = 'pKa',
    str_save = f'reweight_bias.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (-10, 70),
    tup_colormap = (290, 390)
)
#'''
#'''
str_tmp = 'dist_vp_c_pka'
grid_plt(
    list2d_file = [
        [f'290K.0/reweight_bias/{str_tmp}.csv', '290K'],
        [f'310K.0/reweight_bias/{str_tmp}.csv', '310K'],
        [f'330K.0/reweight_bias/{str_tmp}.csv', '330K'],
        [f'350K.0/reweight_bias/{str_tmp}.csv', '350K'],
        [f'370K.0/reweight_bias/{str_tmp}.csv', '370K'],
        [f'390K.0/reweight_bias/{str_tmp}.csv', '390K'],
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
    str_save = f'reweight_bias.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (-10, 70),
    tup_colormap = (290, 390)
)
#'''
'''
str_tmp = 'dist_vp_c_fes'
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
    str_save = f'reweight_bias.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (-10, 70),
    tup_colormap = (290, 390)
)
#'''
plt.show()


