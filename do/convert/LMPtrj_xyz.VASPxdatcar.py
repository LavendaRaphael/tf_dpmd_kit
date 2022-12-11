import ase.io
import numpy

ase_atoms = ase.io.read(
    filename = 'traj.lammpstrj',
    format = 'lammps-dump-text',
    index = ':',
    specorder = ['O','H','C']
)
print(ase_atoms[0])

atoms_xyz = ase.io.read( 
    filename = 'plm.xyz',
    format = 'xyz',
    index = ':'
)
print(atoms_xyz[0])

assert(len(ase_atoms)==len(atoms_xyz))

for int_i in range(len(ase_atoms)):
    atoms_xyz[int_i].set_chemical_symbols('He')
    ase_atoms[int_i].extend(atoms_xyz[int_i])
print(ase_atoms[0])

ase.io.write( filename='XDATCAR', images=ase_atoms, format='vasp-xdatcar' )


