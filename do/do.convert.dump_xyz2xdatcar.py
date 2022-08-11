import ase.io
import numpy
import dpdata
from ase import Atoms
from ase.atom import Atom

type_map = ["O", "H", "C"]
str_origin = "traj.lammpstrj"
dp_sys = dpdata.System(str_origin, fmt='lammps/dump', type_map=type_map)

print(dp_sys)
atoms = dp_sys.to('ase/structure')

atoms_xyz = ase.io.read( 
    filename = 'plm.xyz',
    format = 'xyz',
    index = ':'
)

for int_i in range(len(atoms)):
    atom_vp = Atom(
        symbol = 'He',
        position = atoms_xyz[int_i].get_positions()[6]
    )
    atom_vh1 = Atom(
        symbol = 'Li',
        position = atoms_xyz[int_i].get_positions()[7]
    )
    atom_vh2 = Atom(
        symbol = 'Li',
        position = atoms_xyz[int_i].get_positions()[8]
    )
    atoms[int_i].append(atom_vp)
    atoms[int_i].append(atom_vh1)
    atoms[int_i].append(atom_vh2)

ase.io.write( filename='XDATCAR', images=atoms, format='vasp-xdatcar' )
