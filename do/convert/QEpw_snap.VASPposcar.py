import ase.io
import glob
import os

#'''
str_prefix = 'task.'
str_pwin = 'input'
#'''
'''
str_prefix = 'snap_'
str_pwin = 'pwscf.in'
#'''

list_in = glob.glob( str_prefix + '*/' + str_pwin )

for str_in in list_in:
    str_dir = os.path.dirname(str_in)
    str_poscar = str_dir+'.POSCAR'
    if os.path.isfile( str_poscar ):
        continue
    else:
        ase_atoms = ase.io.read(
            filename = str_in,
            format = 'espresso-in',
        )
        ase.io.write(
            filename = str_poscar,
            images = ase_atoms,
            format = 'vasp',
        )

