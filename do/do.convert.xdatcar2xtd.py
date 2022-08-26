import ase.io

ase_atoms = ase.io.read(
    filename = 'XDATCAR',
    format = 'vasp-xdatcar',
    index=':'
)
ase.io.write( filename='traj.xtd', images=ase_atoms, format='xtd' )
