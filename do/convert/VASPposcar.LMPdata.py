import ase.io

ase_atoms = ase.io.read(filename='POSCAR', format='vasp')
print(ase_atoms)
ase.io.write(filename='lmp.data', format='lammps-data', images=ase_atoms, specorder=['O', 'H', 'C'])
