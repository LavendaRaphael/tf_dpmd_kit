import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def def_plt(
    str_x: list,
    str_y: list,
    list2d_range = None,
    str_save = None
) -> None:
    """

    """
    fig, axe = plt.subplots(1, 1)

    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
    #print(data.dtype)

    h, xedges, yedges, image = axe.hist2d(
        data[str_x],
        data[str_y],
        bins = 100,
        range = list2d_range,
        density = True,
        cmin = 0.001,
        norm = LogNorm()
    )
    fig.colorbar( mappable=image )
    axe.set_xlabel(dict_label[str_x])
    axe.set_ylabel(dict_label[str_y])
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

dict_label = {
    'dist_vp_o_1': 'R(O\N{SUBSCRIPT ONE}V\N{Latin Subscript Small Letter P})',
    'dist_vp_o_2': 'R(O\N{SUBSCRIPT TWO}V\N{Latin Subscript Small Letter P})',
    'del_dist_vp_o': '\N{Greek Capital Letter Delta}R(OV\N{Latin Subscript Small Letter P})',
    'dist_vp_c': 'R(CV\N{Latin Subscript Small Letter P})',

    'distx_vp_o_1': None,
    'distx_vp_o_2': None,
    'cnx_vp_o_1': None,
    'cnx_vp_o_2': None,

    'cn_o_1_h': None,
    'cn_o_2_h': None,
    'del_cn_o_h': '\N{Greek Capital Letter Delta}cn(OH)',

    'cost_o_1_vh': None, 
    'cost_o_2_vh': None,
    'cost_o_h': 'cos dh(OH)',
    'cost_o_1_h': None,
    'cost_o_2_h': None,

    'dist_o_0_h': 'R(O\N{SUBSCRIPT ZERO}H\N{SUBSCRIPT ZERO})',
    'cost_o_0_h': 'cos dh(O\N{SUBSCRIPT ZERO}H)', 

    'pathspath': None,
    'pathzpath': None,
}
for str_key in dict_label.keys():
    if not dict_label[str_key]:
        dict_label[str_key] = str_key

'''
def_plt(
    str_x = 'dist_vp_o_1',
    str_y = 'cn_o_1_h',
    #list2d_range = None,
    #str_save = 'cv.prob.2d.dist_vp_c.del_dist.pdf'
)
#'''
'''
def_plt(
    str_x = 'dist_vp_o_2',
    str_y = 'cn_o_2_h',
    #list2d_range = None,
    #str_save = 'cv.prob.2d.dist_vp_c.del_dist.pdf'
)
#'''
#'''
def_plt(
    str_x = 'dist_vp_c',
    str_y = 'del_dist_vp_o',
    #list2d_range = None,
    str_save = 'cv.prob2d.dist_vp_c.del_dist.pdf'
)
#'''
'''
def_plt(
    str_x = 'dist_vp_o_1',
    str_y = 'dist_vp_o_2',
    #list2d_range = None,
    str_save = 'cv.prob.dist_vp_o_2d.pdf'
)
#'''
plt.show()
