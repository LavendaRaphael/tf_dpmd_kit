import dpdata

dp_sys = dpdata.System('traj.lammpstrj', fmt='lammps/dump', type_map=['O','H','C'])
print(dp_sys)
dp_sys.to('vasp/poscar', 'POSCAR')
