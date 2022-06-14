import dpdata
import ase.io

dp_sys = dpdata.System(
    file_name = 'cp',
    fmt = 'qe/cp/traj',
    )
print(dp_sys)
print(dp_sys['cells'][0])
ase_atoms = dp_sys.to('ase/structure')
ase.io.write( filename='XDATCAR', images=ase_atoms, format='vasp-xdatcar' )
