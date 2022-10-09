import dpdata

dp_sys = dpdata.LabeledSystem(
    file_name = ['input','output'],
    fmt = 'qe/pw/scf'
)
print(dp_sys)

dp_sys.to( 'deepmd/npy', 'dpnpy')


