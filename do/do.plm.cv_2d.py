import numpy as np
import matplotlib.pyplot as plt

def def_plt(
    list2d_cv: list[list],
) -> None:
    """

    """
    int_nplot = len(list2d_cv)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')
    if int_nplot==1:
        axs = [axs]
    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
    # print(data.dtype)
    for int_i in range(int_nplot):
        axe = axs[int_i]
        cvx = list2d_cv[int_i][0]
        cvy = list2d_cv[int_i][1]
        axe.scatter(
            x = data[cvx],
            y = data[cvy]
        )
        axe.set_xlabel( dict_label[cvx] )
        axe.set_ylabel( dict_label[cvy] )
        str_save = list2d_cv[int_i][2]
        if str_save:
            fig.savefig(str_save, bbox_inches='tight')

dict_label = {
    'dist_vp_o_1': None,
    'dist_vp_o_2': None,
    'dist_vp_c': None,

    'distx_vp_o_1': 'R\N{Asterisk}(O\N{SUBSCRIPT ONE}V\N{Latin Subscript Small Letter P})',
    'distx_vp_o_2': 'R\N{Asterisk}(O\N{SUBSCRIPT TWO}V\N{Latin Subscript Small Letter P})',

    'cn_o_1_h': None,
    'cn_o_2_h': None,
    'del_cn_o_h': None,

    'cost_o_1_vh': None, 
    'cost_o_2_vh': None,
    'cost_o_h': 'cos dh(OH)',
    'cost_o_1_h': None,
    'cost_o_2_h': None,

    'dist_o_0_h': 'R(O\N{SUBSCRIPT ZERO}H)',
    'cost_o_0_h': 'cos dh(O\N{SUBSCRIPT ZERO}H)', 

    'pathspath': None,
    'pathzpath': None
}
for str_key in dict_label.keys():
    if not dict_label[str_key]:
        dict_label[str_key] = str_key
#'''
def_plt(
    list2d_cv = [
        #['distx_vp_o_1','distx_vp_o_2','distx.pdf'],
        #['distx_vp_o_1','cost_o_h','distx1_cost.pdf'],
        ['distx_vp_o_2','cost_o_h','distx2_cost.pdf'],
        #['pathspath','pathzpath','path.pdf'],
    ]
)
#'''
plt.show()
