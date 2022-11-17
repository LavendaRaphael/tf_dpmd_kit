import ase.io
import numpy as np
from tf_dpmd_kit import convert

list_ase = ase.io.read(
    filename = 'traj.lammpstrj',
    format = 'lammps-dump-text',
    index = ':',
    specorder = ['O','H','C']
)
print(list_ase[0])

np_snap = range(0, 200000, 500)

convert.list_ase2pwscf(
    list_ase = list_ase,
    np_snap = np_snap
)
