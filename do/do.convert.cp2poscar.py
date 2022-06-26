import dpdata

dp_sys = dpdata.System(
    file_name = 'cp',
    fmt = 'qe/cp/traj',
    )
print(dp_sys)
dp_sys.to(
    fmt = 'vasp/poscar',
    file_name = 'POSCAR',
    frame_idx = -1
    )
