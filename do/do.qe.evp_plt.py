import numpy
from matplotlib import pyplot
import numpy
import matplotlib

np_evp = numpy.genfromtxt('cp/cp.evp')
#np_evp = numpy.genfromtxt('nacl.evp')

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

fig, ax = pyplot.subplots()

ax.plot( np_evp[23:,1], np_evp[23:,2], label = 'Ekin_e', linewidth=2)
ax.plot( np_evp[23:,1], (np_evp[23:,7]-np_evp[23:,5])/10, label = 'Ekin_ion/10', linewidth=2)

ax.legend()
ax.set_xlabel('Time (ps)')
ax.set_ylabel(f'Ekin ($E_h$)')
#ax.set_xlim((0,600))
#ax.set_ylim((0,700))
fig.set_size_inches(8, 4)
fig.savefig('t_ekin.pdf', bbox_inches='tight')
pyplot.show()

