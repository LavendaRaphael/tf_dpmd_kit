from ase import io

ase_atoms = io.read('H2CO3_TS.xsd')
print(ase_atoms)
io.write('H2CO3_TS.pdb', ase_atoms, format='proteindatabank')
