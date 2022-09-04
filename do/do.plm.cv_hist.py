import numpy as np
import matplotlib.pyplot as plt

def def_plt(
    list2d_header: list[list],
    str_save = None
) -> None:
    """

    """
    int_nplot = len(list2d_header)
    fig, axs = plt.subplots(int_nplot, 1)
    if int_nplot==1:
        axs = [axs]

    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
    #print(data.dtype)

    for int_i in range(int_nplot):
        str_header = list2d_header[int_i][0]
        tup_range = list2d_header[int_i][1]
        axs[int_i].hist(
            data[str_header], 
            label = dict_label[str_header],
            bins = 100,
            range = tup_range,
            density = True
        )
        axs[int_i].legend()
    axs[-1].set_xlabel('CV')
    axs[1].set_ylabel('Probability Density')
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

    'cnn_o_w_hmin': r'cn(O$_W$H)$_{min}$'
}

#'''
def_plt(
    list2d_header = [
        ['dist_o_0_h', (0.9, 1.3)],
        ['cost_o_0_h', (-1, -0.9)],
        ['cnn_o_w_hmin', (1.8,1.9)]
    ],
    str_save = 'cv.prob.wall.pdf'
)
#'''
'''
def_plt(
    list2d_header = [
        ['dist_vp_c', (1.2,1.5)],
        ['del_dist_vp_o', (2.1,2.4)],
        ['cost_o_h', (0.8,1)],
    ],
    str_save = 'cv.prob.dist_del_cost.pdf'
)
#'''
plt.show()
