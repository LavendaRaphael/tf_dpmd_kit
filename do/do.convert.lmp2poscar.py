import dpdata

type_map = ["O", "H"]

dp_sys = dpdata.System('water.lmp', fmt='lammps/lmp', type_map=type_map)
print(dp_sys)

dp_sys.to('vasp/poscar', 'POSCAR')
