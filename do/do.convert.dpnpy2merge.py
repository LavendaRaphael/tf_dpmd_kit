import dpdata
import glob

dir_name = './'
file_name = 'type.raw'
target_file_list = sorted(glob.glob('./{}/**/{}'.format(dir_name, file_name), recursive=True))
dp_multisys = dpdata.MultiSystems()
for target_file in target_file_list:
    dp_sys = dpdata.LabeledSystem( target_file[:-8], fmt='deepmd/npy' )
    #dp_sys = dpdata.LabeledSystem( target_file[:-8], fmt='deepmd/raw' )
    print(target_file[:-8])
    print(dp_sys)
    #print(dp_sys['energies'])
    dp_multisys.append(dp_sys)
print(dp_multisys)

dp_multisys.to_deepmd_npy("./")
