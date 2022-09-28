import ase.io
import numpy

ase_atoms = ase.io.read(
    filename='XDATCAR',
    format='vasp-xdatcar',
    index=':'
)
list_cell = ase_atoms[0].get_cell().cellpar()
ase.io.write(
    filename='traj.xyz', 
    images=ase_atoms, 
    format='extxyz', 
    comment=f' {list_cell[0]} {list_cell[1]} {list_cell[2]}' 
)
