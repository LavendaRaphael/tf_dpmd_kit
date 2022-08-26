import ase.io

ase_atoms = ase.io.read(
    filename = 'XDATCAR',
    format = 'vasp-xdatcar',
    index=':'
)
print(ase_atoms[0].get_cell())
print(ase_atoms[1].get_cell())
ase.io.write( filename='traj.pdb', images=ase_atoms, format='proteindatabank' )
