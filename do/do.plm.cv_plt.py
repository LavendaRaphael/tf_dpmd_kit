import numpy as np
import matplotlib.pyplot as plt

def def_plt(
    list2d_label: list[list],
) -> None:
    """

    """
    int_nplot = len(list2d_label)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')
    
    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
    #print(data.dtype)
    
    for int_i in range(int_nplot):
        list_label = list2d_label[int_i]
        for str_header in list_label:
            axs[int_i].plot(data['time'], data[str_header], label=dict_label[str_header])
            axs[int_i].legend() 
    axs[-1].set_xlabel('Time(ps)')
    axs[0].set_ylabel('CV')

dict_label = {
    "dist_vp_c": 'R(CP)',
    "cn_o_1_h": None,
    'cn_o_2_h': None,
    'del_cn_o_h': 'cn(O\N{SUBSCRIPT ONE}H)-cn(O\N{SUBSCRIPT TWO}H)',
    'cn_o_1_vp': None,
    'cn_o_2_vp': None,
    'cost_o_1_vp': None, 
    'cost_o_2_vp': None,
    'cost_o_1_vh': None, 
    'cost_o_2_vh': None, 
    'cost_o_1_h': 'cos dh(O\N{SUBSCRIPT ONE}H)', 
    'cost_o_2_h': 'cos dh(O\N{SUBSCRIPT TWO}H)', 
    'restraintbias': None,
    'up_c_vpbias': None,
    'up_c_vp.force2': None, 
    'dist_o_0_h': 'R(O\N{SUBSCRIPT ZERO}H)',
    'up_o_0_hbias': None, 
    'up_o_0_h.force2': None, 
    'cost_o_0_h': 'cos dh(O\N{SUBSCRIPT ZERO}H)', 
    'dist_vp_o_1': None, 
    'dist_vp_o_2': None, 
    'dist_o_1_vh': None, 
    'dist_o_2_vh': None
}

def_plt(
    list2d_label = [
        ["dist_vp_c"],
        ['del_cn_o_h'],
        ['cost_o_1_h','cost_o_2_h','cost_o_0_h'],
        ['dist_o_0_h'],
    ]
)
plt.show()
