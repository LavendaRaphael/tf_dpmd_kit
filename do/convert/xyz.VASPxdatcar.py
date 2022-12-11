import ase.io
import numpy as np

list_atoms = ase.io.read(
    filename = 'plm.xyz',
    format = 'xyz',
    index = ':'
)
print(list_atoms[0])

for atoms in list_atoms:
    atoms.set_cell( 15.679309*np.identity(3) )
    atoms.set_pbc((True, True, True))

ase.io.write( filename='XDATCAR', images=list_atoms, format='vasp-xdatcar' )
