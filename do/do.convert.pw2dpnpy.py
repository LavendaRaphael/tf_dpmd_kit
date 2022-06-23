import dpdata

dp_sys = dpdata.MultiSystems.from_dir(
    dir_name = './',
    file_name = 'pwscf.out',
    fmt = 'qe/pw/scf'
    )
print(dp_sys)

dp_sys.to( 'deepmd/npy', 'dpnpy')


