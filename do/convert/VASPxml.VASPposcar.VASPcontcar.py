import ase.io

list_atoms = ase.io.read(
    filename = 'vasprun.xml',
    index = ':',
    format = 'vasp-xml',
)
print(list_atoms)
ase.io.write(
    filename = 'POSCAR',
    images = list_atoms[0],
    format = 'vasp',
)
ase.io.write(
    filename = 'CONTCAR',
    images = list_atoms[-1],
    format = 'vasp',
)
