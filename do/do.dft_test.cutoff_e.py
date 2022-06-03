import numpy
import os
import ase.io
from matplotlib import pyplot
import matplotlib

'''
np_cutoff = numpy.arange(
    start = 700,
    stop = 1600,
    step = 100,
    )
str_filename = 'OUTCAR'
str_format = 'vasp-out'
'''
np_cutoff = numpy.arange(
    start = 100,
    stop = 190,
    step = 10,
    )
str_filename = 'log'
str_format = 'espresso-out'
#'''

np_energy = numpy.zeros( shape = np_cutoff.shape )
for int_i in range( np_cutoff.shape[0] ):
    int_cutoff = np_cutoff[ int_i ]
    print( int_cutoff )
    os.chdir( str(int_cutoff) )
    np_energy[ int_i ] = ase.io.read(
        filename = str_filename,
        format = str_format
        ).get_total_energy()
    os.chdir('..')

print(np_energy)
np_energy -= np_energy[-1]
print(np_energy)

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

fig, ax = pyplot.subplots()

ax.plot( np_cutoff, np_energy, label = '', linewidth=2)

ax.legend()
ax.set_xlabel('')
ax.set_ylabel('')
#ax.set_xlim((0,600))
#ax.set_ylim((0,700))
fig.set_size_inches(7, 6)
fig.savefig('converge_cutoff.pdf', bbox_inches='tight')
pyplot.show()
