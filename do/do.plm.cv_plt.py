import numpy as np
import matplotlib.pyplot as plt

def def_plt(
    list2d_header: list[list],
    dict_windows: dict
) -> None:
    """

    """
    int_nplot = len(list2d_header)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')
    if int_nplot==1:
        axs = [axs]

    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header, invalid_raise=False)
    #print(data.dtype)

    data_windows = data[ dict_windows['list_range'][0]:dict_windows['list_range'][1] ]

    for int_i in range(int_nplot):
        list_header = list2d_header[int_i]
        for str_header in list_header:
            if str_header not in dict_label:
                str_label = str_header
            else:
                str_label = dict_label[str_header]
            axs[int_i].plot(data_windows['time']/dict_windows['timestep'], data_windows[str_header], label=str_label)
            axs[int_i].legend() 
    #axs[-1].set_xlabel('Time(ps)')
    axs[0].set_ylabel('CV')

dict_label = {
    'dist_vp_o_1': 'R(O\N{SUBSCRIPT ONE}V\N{Latin Subscript Small Letter P})',
    'dist_vp_o_2': 'R(O\N{SUBSCRIPT TWO}V\N{Latin Subscript Small Letter P})',
    'del_dist_vp_o': '\N{Greek Capital Letter Delta}R(OV\N{Latin Subscript Small Letter P})',
    'dist_vp_c': 'R(CV\N{Latin Subscript Small Letter P})',
    'cost_o_h': 'cos dh(OH)',
    'dist_o_0_h': 'R(O\N{SUBSCRIPT ZERO}H)',
    'cost_o_0_h': 'cos dh(O\N{SUBSCRIPT ZERO}H)', 
}
#'''
def_plt(
    list2d_header = [
        ['dist_vp_c'],
        ['cost_o_0_h'],
        ['metadbias'],
        ['wall_vp_cbias'],
        ['dist_o_0_h'],
        ['wall_dist_o_0_hbias'],
        ['cn_o_0_h'],
        ['wall_cn_o_0_hbias'],
        ['cnn_o_w_hmin'],
        ['wall_cnn_o_w_hbias']
    ],
    dict_windows = {
        'list_range': [0, 2000000],
        'timestep': 0.0005
    }
)
plt.show()
