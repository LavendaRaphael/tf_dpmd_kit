import numpy
import os
import ase.io
from matplotlib import pyplot
import matplotlib
from numpy import linalg as LA
import dpdata

#'''
list_std = ['2000/OUTCAR', 'vasp-out']
list_file = ['1900/OUTCAR', 'vasp-out', 'VASP(1900eV) - VASP(2000eV)']
str_filesave = 'atom_1900.pdf'
#'''
'''
list_std = ['300/log', 'espresso-out']
list_file = ['280/log', 'espresso-out', 'QE(280Ry) - QE(300Ry)']
str_filesave = 'atom_280.pdf'
#'''

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

np_force_std = ase.io.read(
    filename = list_std[0],
    format = list_std[1]
    ).get_forces()

np_atoms = numpy.arange( np_force_std.shape[0] )

fig, ax = pyplot.subplots()
np_force = ase.io.read(
    filename = list_file[0],
    format = list_file[1]
    ).get_forces()

np_force_error = np_force - np_force_std

np_force_error_norm = LA.norm(np_force_error, axis=1)
print(np_force_error_norm)

ax.plot( np_atoms, np_force_error_norm, label = list_file[2], linewidth=2)

ax.legend()
ax.set_xlabel( 'Atom' )
ax.set_ylabel( 'Force error (eV/Å)' )
#ax.set_xlim((0,600))
#ax.set_ylim((0,700))
fig.set_size_inches(7, 6)
fig.savefig(str_filesave, bbox_inches='tight')
pyplot.show()
