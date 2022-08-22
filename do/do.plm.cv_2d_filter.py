import numpy as np
import matplotlib.pyplot as plt

def def_plt(
    list2d_cv: list[list],
    list2d_filter: list[list]
) -> None:
    """

    """
    int_nplot = len(list2d_cv)
    fig, axs = plt.subplots(int_nplot, 1)
    if int_nplot==1:
        axs = [axs]

    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
    # print(data.dtype)

    np_filter = np.empty( shape=(data.shape[0]) )
    np_filter = True
    for list_filter in list2d_filter:
        str_label = list_filter[0]
        float_min = list_filter[1]
        float_max = list_filter[2]
        np_filter = np_filter & (data[str_label] > float_min) & (data[str_label] < float_max)
    data_filter = data[np_filter]

    for int_i in range(int_nplot):
        axe = axs[int_i]
        list_cv = list2d_cv[int_i]
        cvx = list_cv[0]
        cvy = list_cv[1]
        str_save = list_cv[2]
        axe.scatter(
            x = data[cvx],
            y = data[cvy]
        )
        axe.scatter(
            x = data_filter[cvx],
            y = data_filter[cvy]
        )
        axe.set_xlabel( dict_label[cvx] )
        axe.set_ylabel( dict_label[cvy] )
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
        ['distx_vp_o_1','distx_vp_o_2', None],
        ['distx_vp_o_1','cost_o_h', None],
        ['distx_vp_o_2','cost_o_h', None],
        ['pathspath','pathzpath', None],
    ],
    #list2d_filter = [
    #    ['distx_vp_o_1', 0, 1],
    #    ['distx_vp_o_2', 0, 4],
    #    ['cost_o_h', 0.75, 1.5]
    #]
    #list2d_filter = [
    #    ['distx_vp_o_1', 0, 4],
    #    ['distx_vp_o_2', 0, 1],
    #    ['cost_o_h', 0.75, 1.5]
    #]
    # C+C
    #list2d_filter = [
    #    ['distx_vp_o_1', 1.5, 3],
    #    ['distx_vp_o_2', 3.5, 5],
    #    ['cost_o_h', -1, 2]
    #]
    # bulk
    list2d_filter = [
        ['distx_vp_o_1', 3.8, 5],
        ['distx_vp_o_2', 3.8, 5],
        ['cost_o_h', -1, 2]
    ]
)
#'''
plt.show()
