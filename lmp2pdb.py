import ase.io
import dpdata

dp_sys = dpdata.System()
for int_i in range(0,10000,40):
    str_origin = "{}.lammpstrj".format(int_i)
    str_new = "{:04d}".format(int_i)
    dp_tmp = dpdata.System(str_origin, fmt='lammps/dump', type_map=["O","H"])
    dp_sys.append(dp_tmp)
print(dp_sys)
ase_atoms = dp_sys.to('ase/structure')
ase.io.write(format="proteindatabank", filename="traj.pdb", images=ase_atoms)
