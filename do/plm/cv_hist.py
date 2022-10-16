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
        ax = axs[int_i]
        if str_header not in dict_label:
            str_label = str_header
        else:
            str_label = dict_label[str_header]
        ax.hist(
            data[str_header],
            label = str_label,
            bins = 'auto',
            density = True
        )
        if len(list2d_header[int_i]) > 1:
            ax.set_xlim(list2d_header[int_i][1])

        float_std = np.std(data[str_header])
        float_mean = np.mean(data[str_header])
        ax.plot([],[],' ',label=f'mean = {float_mean:.3f}')
        ax.plot([],[],' ',label=f'std = {float_std:.3f}')

        ax.legend()
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
        ['dist_vp_o_1_2', (1.05,1.25)],
        ['dist_o_0_h', (0.9, 1.2)],
        ['cn_o_0_h', (0,0.01)],
    ],
    str_save = 'cv.prob.pdf'
)
#'''
plt.show()
