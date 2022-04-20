from ase import Atoms
from deepmd.calculator import DP
import dpdata
from numpy import linalg as LA
from matplotlib import pyplot as plt 
import numpy as np

list_model = [
    "../../../00.train/000/frozen_model.pb",
    "../../../00.train/001/frozen_model.pb",
    "../../../00.train/002/frozen_model.pb",
    "../../../00.train/003/frozen_model.pb" 
    ]

dp_sys = dpdata.System("../traj/3800.lammpstrj", fmt='lammps/dump', type_map=["O","H"])
ase_atoms = dp_sys.to('ase/structure')[0]
print(ase_atoms)

int_natom = len(ase_atoms.get_atomic_numbers())
print(int_natom)

array3d_forces = np.zeros(shape=(4,int_natom,3))
for int_i in range(4):
    ase_atoms.calc = DP( model=list_model[int_i] )
    array3d_forces[int_i] = ase_atoms.get_forces()

array2d_forces_ave = np.average (array3d_forces, axis=0)

array3d_devi_f = np.zeros(shape=(4,int_natom,3))
for int_i in range(4):
    array3d_devi_f[int_i] = array3d_forces[int_i] - array2d_forces_ave

array2d_devi_f_norm = LA.norm(array3d_devi_f, axis=2)**2

array1d_devi_f_final = np.sqrt(np.average (array2d_devi_f_norm, axis=0))


float_lowthred = 0.15
for int_i in range(int_natom):
    if array1d_devi_f_final[int_i] > float_lowthred:
        print(int_i, array1d_devi_f_final[int_i])

x = np.arange(int_natom)
fig, ax = plt.subplots()
ax.plot(x, array1d_devi_f_final)
ax.legend()
plt.show()
