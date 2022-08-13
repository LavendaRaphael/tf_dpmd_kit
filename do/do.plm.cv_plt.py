import numpy as np
import matplotlib.pyplot as plt

def def_plt(
    list_label: list,
    str_save: str = None
) -> None:
    """

    """
    fig, ax = plt.subplots()
    
    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
    #print(data.dtype)
    
    for str_header in list_label:
        ax.plot(data['time'], data[str_header], label=dict_label[str_header])
    
    ax.legend()
    ax.set_xlabel('Time(ps)')
    ax.set_ylabel('CV')
    if str_save:
        ax.savefig(str_save, bbox_inches='tight')
    #plt.show()

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
    list_label = [
        "dist_vp_c",
    ],
    #str_save = 'cv.dist_vp_c.pdf'
)

def_plt(
    list_label = [
        'del_cn_o_h',
    ],
    #str_save = 'cv.del_cn_o_h.pdf'
)
def_plt(
    list_label = [
        'cost_o_1_h', 
        'cost_o_2_h', 
        'cost_o_0_h', 
    ],
    #str_save = 'cv.cost.pdf'
)
def_plt(
    list_label = [
        'dist_o_0_h',
    ],
    #str_save = 'cv.dist_o_0_h.pdf'
)
plt.show()
