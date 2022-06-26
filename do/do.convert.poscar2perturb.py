import dpdata
import os

dp_sys = dpdata.System(
    file_name = 'POSCAR',
    fmt = 'vasp/poscar',
    )

int_perturb = 100
dp_sys_perturb = dp_sys.perturb(
    pert_num = int_perturb,
    cell_pert_fraction = 0.0,
    atom_pert_distance = 0.03,
    )

str_dir = 'perturb'
if not os.path.exists(str_dir):
    os.mkdir(str_dir)

for int_i in range(int_perturb):
    dp_sys_perturb.to( 'vasp/poscar', os.path.join(str_dir, f'{int_i:02d}.POSCAR'), frame_idx=int_i )
