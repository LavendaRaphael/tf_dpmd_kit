import numpy
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import rc
import os
from matplotlib import cm

def grid_plt(
    list2d_file: list[list],
    str_xlabel: str,
    str_save: str=None,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    tup_colormap: tuple = None
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
        ax.plot( np_data[:,0], np_data[:,1]-min(np_data[:,1]), label=str_label, linewidth=2, color=color)
    
    ax.legend()
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel('FES (kJ/mol)')
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.set_size_inches(9, 7)
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
#'''
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
    tup_ylim = (-10, 70),
    tup_colormap = (290, 390)
)
#'''
plt.show()


