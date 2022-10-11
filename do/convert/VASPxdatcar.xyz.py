import ase.io
import numpy

ase_atoms = ase.io.read(
    filename='XDATCAR',
    format='vasp-xdatcar',
    index=':'
)
ase.io.write(
    filename='traj.xyz', 
    images=ase_atoms, 
    format='xyz' 
)
