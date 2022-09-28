import ase.io
import numpy
from ase import Atoms
from ase.atom import Atom
from tf_dpmd_kit import convert

array_id = numpy.arange(0,5005,5)

atoms = convert.def_dump2ase(
    array_id = array_id,
    type_map = ["O", "H", "C"]
)

atoms_xyz = ase.io.read( 
    filename = 'plm.xyz',
    format = 'xyz',
    index = ':'
)

for int_i in range(len(atoms)):
    atom_vp = Atom(
        symbol = 'He',
        position = atoms_xyz[int_i].get_positions()[0]
    )
    atom_vh1 = Atom(
        symbol = 'Li',
        position = atoms_xyz[int_i].get_positions()[1]
    )
    atom_vh2 = Atom(
        symbol = 'Li',
        position = atoms_xyz[int_i].get_positions()[2]
    )
    atoms[int_i].append(atom_vp)
    atoms[int_i].append(atom_vh1)
    atoms[int_i].append(atom_vh2)

ase.io.write( filename='XDATCAR', images=atoms, format='vasp-xdatcar' )
