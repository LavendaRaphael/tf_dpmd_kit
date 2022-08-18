import dpdata
import ase.io
import numpy

dp_sys = dpdata.System('traj.lammpstrj', fmt='lammps/dump', type_map=['O','H','C'])
print(dp_sys)
ase_atoms = dp_sys.to('ase/structure')
list_cell = ase_atoms[0].get_cell().cellpar()
ase.io.write( filename='traj.xyz', images=ase_atoms, format='extxyz', comment=f' {list_cell[0]} {list_cell[1]} {list_cell[2]}' )
