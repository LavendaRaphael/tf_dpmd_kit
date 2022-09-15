import numpy
from matplotlib import pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

def grid_plt(
    list2d_file: list[list],
    str_save: str=None
):
    fig, ax = plt.subplots()
    for list_file in list2d_file:
        str_file = list_file[0]
        str_label = list_file[1]

        with open(str_file, 'r') as open_file:
            list_header = open_file.readline().split()[2:]
        np_data = np.genfromtxt(str_file, names=list_header)
        
        ax.plot( np_data['dist_vp_c'], np_data['fes_dist_vp_c'], label=str_label)
    
    ax.legend()
    ax.set_xlabel('R(CP) (Å)')
    ax.set_ylabel('FES (kJ/mol)')
    #ax.set_ylim()
    #fig.set_size_inches(8, 4)
    if str_save:
        fig.savefig('fes.dist_vp_c.pdf', bbox_inches='tight')

grid_plt(
    list2d_file = [
        ['analysis.1.fes.dist_vp_c.grid', '0.2ns'],
        ['analysis.3.fes.dist_vp_c.grid', '0.4ns'],
        ['analysis.4.fes.dist_vp_c.grid', '0.6ns'],
        ['fes.dist_vp_c.grid', '0.8ns'],
    ],
    str_save = 'fes.dist_vp_c.pdf'
)

plt.show()


