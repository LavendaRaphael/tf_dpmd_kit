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
    filename = 'plm.vp.xyz',
    format = 'xyz',
    index = ':'
)
print(atoms_xyz[0])

assert(len(ase_atoms)==len(atoms_xyz))

for int_i in range(len(ase_atoms)):
    atoms_xyz[int_i].set_chemical_symbols('He')
    ase_atoms[int_i].extend(atoms_xyz[int_i])
print(ase_atoms[0])

list_cell = ase_atoms[0].get_cell().cellpar()
ase.io.write( filename='traj.xyz', images=ase_atoms, format='extxyz', comment=f' {list_cell[0]} {list_cell[1]} {list_cell[2]}' )
