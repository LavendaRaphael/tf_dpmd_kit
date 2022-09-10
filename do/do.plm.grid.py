import numpy
from matplotlib import pyplot
import numpy as np
import matplotlib

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

str_evp = 'fes_dist_vp_c.grid'
with open(str_evp, 'r') as open_evp:
    list_header = open_evp.readline().split()[2:]
np_evp = np.genfromtxt(str_evp, names=list_header)
print(np_evp.dtype)

fig, ax = pyplot.subplots()
ax.plot( np_evp['dist_vp_c'], np_evp['fes_dist_vp_c'], linewidth=2)

#ax.legend()
ax.set_xlabel('R(CP) (Å)')
ax.set_ylabel('FES (kJ)')
#ax.set_ylim()
#fig.set_size_inches(8, 4)
fig.savefig('fes.dist_vp_c.pdf', bbox_inches='tight')
pyplot.show()

