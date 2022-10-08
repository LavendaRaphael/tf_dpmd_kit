import ase.io

list_ase = ase.io.read(
    filename = 'traj.lammpstrj',
    format = 'lammps-dump-text',
    index = ':',
    specorder = ['O','H','C']
)
print(list_ase[0])

ase.io.write( filename='traj.xyz', images=list_ase, format='xyz' )
