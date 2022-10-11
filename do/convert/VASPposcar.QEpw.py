import ase.io
from tf_dpmd_kit import convert

ase_atoms = ase.io.read(
    filename = 'POSCAR',
    format = 'vasp',
)

convert.ase2pwscf(
    ase_atoms = ase_atoms
)

