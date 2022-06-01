import dpdata

dp_sys = dpdata.System(
    file_name = 'cp',
    fmt = 'qe/cp/traj',
    )
dp_sys.to(
    fmt = 'vasp/poscar',
    file_name = 'POSCAR',
    frame_idx = -1
    )
