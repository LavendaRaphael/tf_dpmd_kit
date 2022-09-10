import numpy
from matplotlib import pyplot as plt
import numpy as np
import matplotlib
from scipy.interpolate import interp2d

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

str_evp = 'fes_cost_o_h.grid'
with open(str_evp, 'r') as open_evp:
    list_header = open_evp.readline().split()[2:]
np_evp = np.genfromtxt(str_evp, names=list_header)
print(np_evp.dtype)

x_list=np_evp['cost_o_0_h']
y_list=np_evp['del_cost_o_h']
z_list=np_evp['fes_cost_o_h']

N = int(len(z_list)**.5)
z = z_list.reshape(N, N)
plt.imshow(z, extent=(np.amin(x_list), np.amax(x_list), np.amin(y_list), np.amax(y_list)), aspect = 'auto')
plt.colorbar()

#fig, ax = pyplot.subplots()
#ax.plot( np_evp['dist_vp_c'], np_evp['fes_dist_vp_c'], linewidth=2)

#ax.legend()
#ax.set_xlabel('R(CP) (Å)')
#ax.set_ylabel('FES (kJ)')
#ax.set_ylim()
#fig.set_size_inches(8, 4)
#fig.savefig('fes.dist_vp_c.pdf', bbox_inches='tight')
plt.show()

