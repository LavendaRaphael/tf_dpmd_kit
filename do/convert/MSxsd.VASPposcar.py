from ase import io

ase_atoms = io.read('HCO3.xsd')
print(ase_atoms)
io.write('HCO3.POSCAR', ase_atoms, format='vasp')
