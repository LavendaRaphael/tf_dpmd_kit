import numpy
from matplotlib import pyplot
import numpy
import matplotlib

np_evp = numpy.genfromtxt('cp/cp.evp')
print(np_evp)

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

fig, ax = pyplot.subplots()

ax.plot( np_evp[23:,1], np_evp[23:,2], label = '', linewidth=2)

ax.legend()
ax.set_xlabel('Time (ps)')
ax.set_ylabel(f'Ekinc ($E_h$)')
#ax.set_xlim((0,600))
#ax.set_ylim((0,700))
fig.set_size_inches(7, 6)
fig.savefig('t_ekinc.pdf', bbox_inches='tight')
pyplot.show()

