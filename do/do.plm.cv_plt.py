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
    dict_label = {
        #"dist_pro_c",
        #"cn_o_1_h",
        #'cn_o_2_h',
        #'del_cn_o_h',
        #'cn_c_pro',
        #'cn_c_pro_x',
        #'cost_o_1_pro', 
        #'cost_o_2_pro',
        #'cost_o_1_vh', 
        #'cost_o_2_vh', 
        #'cost_o_1_h', 
        #'cost_o_2_h', 
        #'restraint.bias',
        #'up_c_pro.bias',
        #'up_c_pro.force2', 
        #'dist_o_0_h',
        #'up_o_0_h.bias', 
        #'up_o_0_h.force2', 
        #'cost_o_0_h', 
        #'dist_pro_o_1', 
        #'dist_pro_o_2', 
        #'dist_o_1_vh', 
        #'dist_o_2_vh'
    }
)
