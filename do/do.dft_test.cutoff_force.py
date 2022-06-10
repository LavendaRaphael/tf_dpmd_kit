import numpy
import os
import ase.io
from matplotlib import pyplot
import matplotlib
from numpy import linalg as LA
import dpdata

#'''
np_cutoff = numpy.arange(
    start = 700,
    stop = 2100,
    step = 100,
    )
str_filename = 'OUTCAR'
str_format = 'vasp-out'
str_xlabel = 'Cutoff (eV)'
list_std = ['2000/OUTCAR', 'vasp-out','converge.pdf', 'VASP - VASP(2000eV)']
#list_std = ['1600/OUTCAR', 'vasp-out','converge.pdf', 'VASP - VASP(1600eV)']
#list_std = ['../qe_cutoff_ppsg15/300/log', 'espresso-out','converge_qe_ppsg15_300.pdf', 'VASP - QE(300Ry)']
#'''
'''
np_cutoff = numpy.array([100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300])
str_filename = 'log'
str_format = 'espresso-out'
str_xlabel = 'Cutoff (Ry)'
list_std = ['300/log', 'espresso-out','converge.pdf', 'QE - QE (300Ry)']
#list_std = ['280/log', 'espresso-out','converge_qe_280.pdf', 'QE - QE (280Ry)']
#list_std = ['../vasp_cutoff/2000/OUTCAR', 'vasp-out','converge_vasp_2000.pdf', 'QE - VASP(2000eV)']
#'''
'''
np_cutoff = numpy.array([100,110,120,130,140,150,160,170,180,190,200])
str_filename = 'log'
str_format = 'espresso-out'
str_xlabel = 'Cutoff (Ry)'
#list_std = ['300/log', 'espresso-out','converge.pdf', 'QE - QE (300Ry)']
list_std = ['../qe_cutoff_ppsg15/300/log', 'espresso-out','converge_qe_ppsg15_300.pdf', 'QE - QE(300Ry)']
#'''


ase_atoms = ase.io.read(
    filename = list_std[0],
    format = list_std[1]
    )
list_symbols = ase_atoms.get_chemical_symbols()
np_filter = (numpy.array(list_symbols) != 'H')

np_force_std = ase_atoms.get_forces()

np_force_error_norm_ave = numpy.zeros( shape = np_cutoff.shape[0] )

for int_i in range( np_cutoff.shape[0] ):
    int_cutoff = np_cutoff[ int_i ]
    print( int_cutoff )
    os.chdir( str(int_cutoff) )

    np_force = ase.io.read(
        filename = str_filename,
        format = str_format
        ).get_forces()

    np_force_error = np_force - np_force_std

    np_force_error_norm = LA.norm(np_force_error, axis=1)

    print(np_force_error_norm[np_filter])
    np_force_error_norm_ave[ int_i ] = numpy.average( np_force_error_norm[np_filter] )

    os.chdir('..')

for int_i in range(np_cutoff.shape[0]):
    print(np_cutoff[int_i], f'{np_force_error_norm_ave[int_i]:.3f}')

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

fig, ax = pyplot.subplots()

ax.plot( np_cutoff, np_force_error_norm_ave, label = list_std[3], linewidth=2)

ax.legend()
ax.set_xlabel( str_xlabel )
ax.set_ylabel( 'Average force devi (eV/Å)' )
#ax.set_xlim((0,600))
ax.set_ylim((0,0.3))
fig.set_size_inches(7, 6)
fig.savefig(list_std[2], bbox_inches='tight')
pyplot.show()
