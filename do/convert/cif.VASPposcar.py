import ase.io

ase_atoms = ase.io.read(
     filename = 'H2CO3_CC.cif',
     format = 'cif',
     )
ase.io.write(
    filename = 'H2CO3_CC.POSCAR',
    images = ase_atoms,
    format = 'vasp',
    )

