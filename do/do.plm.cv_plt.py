import numpy as np
import matplotlib.pyplot as plt

def def_plt(
    list2d_label: list[list],
    dict_windows: dict
) -> None:
    """

    """
    int_nplot = len(list2d_label)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')
    if int_nplot==1:
        axs = [axs]

    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
    #print(data.dtype)

    list_range = [None, None]
    list_range[0] = dict_windows['list_range'][0] // dict_windows['int_stride']
    list_range[1] = dict_windows['list_range'][1] // dict_windows['int_stride']
    data_windows = data[ list_range[0]:list_range[1] ]

    for int_i in range(int_nplot):
        list_label = list2d_label[int_i]
        for str_header in list_label:
            axs[int_i].plot(data_windows['time'], data_windows[str_header], label=dict_label[str_header])
            axs[int_i].legend() 
    axs[-1].set_xlabel('Time(ps)')
    axs[0].set_ylabel('CV')

dict_label = {
    'dist_vp_o_1': None,
    'dist_vp_o_2': None,
    'dist_vp_c': None,

    'distx_vp_o_1': None,
    'distx_vp_o_2': None,
    'cnx_vp_o_1': None,
    'cnx_vp_o_2': None,

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
    'pathzpath': None,
    
    'restraintbias': None
}
for str_key in dict_label.keys():
    if not dict_label[str_key]:
        dict_label[str_key] = str_key
#'''
def_plt(
    list2d_label = [
        ['cn_o_1_h','cn_o_2_h'],
        ['dist_vp_c'],
        ['cost_o_h'],
        #['pathspath'],
        #['pathzpath'],
        ['dist_o_0_h'],
        ['cost_o_0_h'],
        #['cnx_vp_o_1','cnx_vp_o_2'],
        ['restraintbias'],
    ],
    dict_windows = {
        'list_range': [86750,86950],
        'int_stride': 20
    }
)
plt.show()
