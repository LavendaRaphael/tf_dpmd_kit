import numpy as np
import matplotlib.pyplot as plt

def def_plt(
    str_header: str,
    str_save = None
) -> None:
    """

    """
    fig, ax = plt.subplots()

    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
    #print(data.dtype)

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
    #ax.set_xlim(list2d_header[int_i][1])

    float_std = np.std(data[str_header])
    float_mean = np.mean(data[str_header])
    ax.plot([],[],' ',label=f'mean = {float_mean:.3f}')
    ax.plot([],[],' ',label=f'std = {float_std:.3f}')

    ax.legend()
    ax.set_xlabel('CV')
    ax.set_ylabel('Probability Density')
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
list_header = [
    'dist_c_o_0','dist_c_o_1', 'cn_c_o','dist_c_h_0', 'cn_c_h','dist_o_0_h','cn_o_1_2_h'
]
for str_tmp in list_header:
    def_plt(
        str_header = str_tmp,
        str_save = f'{str_tmp}.prob.pdf'
    )
#'''
plt.show()
