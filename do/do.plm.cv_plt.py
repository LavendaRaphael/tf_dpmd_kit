import numpy as np
import matplotlib.pyplot as plt

def def_plt(
    dict_label: list,
    str_save: str = None
) -> None:
    """

    """
    fig, ax = plt.subplots()
    
    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
    #print(data.dtype)
    
    for str_label in dict_label:
        ax.plot(data['time'], data[str_label], label=str_label)
    
    ax.legend()
    ax.set_xlabel('Time(ps)')
    ax.set_ylabel('CV')
    if str_save:
        plt.savefig(str_save, bbox_inches='tight')
    plt.show()

def_plt(
    dict_label = [
        "dist_vp_c",
    ]
)

def_plt(
    dict_label = [
        'del_cn_o_h',
    ]
)
def_plt(
    dict_label = [
        'cost_o_1_h', 
        'cost_o_2_h', 
        'cost_o_0_h', 
    ]
)
def_plt(
    dict_label = [
        'dist_o_0_h',
    ]
)
'''
def_plt(
    dict_label = [
        #"dist_vp_c",
        #"cn_o_1_h",
        #'cn_o_2_h',
        #'del_cn_o_h',
        #'cn_o_1_vp',
        #'cn_o_2_vp',
        #'cost_o_1_vp', 
        #'cost_o_2_vp',
        #'cost_o_1_vh', 
        #'cost_o_2_vh', 
        #'cost_o_1_h', 
        #'cost_o_2_h', 
        #'restraintbias',
        #'up_c_vpbias',
        #'up_c_vp.force2', 
        #'dist_o_0_h',
        #'up_o_0_hbias', 
        #'up_o_0_h.force2', 
        #'cost_o_0_h', 
        #'dist_vp_o_1', 
        #'dist_vp_o_2', 
        #'dist_o_1_vh', 
        #'dist_o_2_vh'
    ]
)
'''
