import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def def_plt(
    str_x: str,
    str_y: str,
    list2d_range = None,
    str_save = None
) -> None:
    """

    """
    fig, ax = plt.subplots(1, 1)

    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header, invalid_raise=False)
    #print(data.dtype)

    h, xedges, yedges, image = ax.hist2d(
        data[str_x],
        data[str_y],
        bins = 100,
        range = list2d_range,
        density = True,
        cmin = 0.001,
        norm = LogNorm()
    )
    fig.colorbar( mappable=image )
    if str_x not in dict_label:
        dict_label[str_x] = str_x
    if str_y not in dict_label:
        dict_label[str_y] = str_y
    ax.set_xlabel(dict_label[str_x])
    ax.set_ylabel(dict_label[str_y])
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

dict_label = {
    'dist_vp_o_1': 'R(O\N{SUBSCRIPT ONE}V\N{Latin Subscript Small Letter P})',
    'dist_vp_o_2': 'R(O\N{SUBSCRIPT TWO}V\N{Latin Subscript Small Letter P})',
    'del_dist_vp_o': '\N{Greek Capital Letter Delta}R(OV\N{Latin Subscript Small Letter P})',
    'dist_vp_c': 'R(CV\N{Latin Subscript Small Letter P})',

    'del_cn_o_h': '\N{Greek Capital Letter Delta}cn(OH)',

    'cost_o_h': 'cos dh(OH)',

    'dist_o_0_h': 'R(O\N{SUBSCRIPT ZERO}H\N{SUBSCRIPT ZERO})',
    'cost_o_0_h': 'cos dh(O\N{SUBSCRIPT ZERO}H)', 
}
for str_key in dict_label.keys():
    if not dict_label[str_key]:
        dict_label[str_key] = str_key

#'''
def_plt(
    str_x = 'dist_vp_c',
    str_y = 'cost_o_0_h',
    #list2d_range = None,
    str_save = 'cv.prob2d.dist_vp_c.cost_o_0_h.pdf'
)
#'''
plt.show()
