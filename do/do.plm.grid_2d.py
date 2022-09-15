import numpy
from matplotlib import pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

def grid2d_plt(
    str_file: str,
    str_save: str=None
):
    fig, ax = plt.subplots()

    with open(str_file, 'r') as open_file:
        list_header = open_file.readline().split()[2:]
    np_data = np.genfromtxt(str_file, names=list_header)
    
    x_list=np_data['dhx_o_0_h']
    y_list=np_data['dhx_o_1_2_vh']
    z_list=np_data['fes_dh_o_h']

    N = int(len(z_list)**.5)
    z = z_list.reshape(N, N)
    image = ax.imshow(z, extent=(np.amin(x_list), np.amax(x_list), np.amin(y_list), np.amax(y_list)), aspect = 'auto')

    fig.colorbar(mappable=image)
    ax.set_xlabel(r'DH(O_0H_0)')
    ax.set_ylabel('DH(OH)')
    #fig.set_size_inches(8, 4)
    if str_save:
        fig.savefig('fes.dist_vp_c.pdf', bbox_inches='tight')

grid2d_plt(
    str_file = 'fes.dh_o_h.grid'
    #str_save = 
)

plt.show()


