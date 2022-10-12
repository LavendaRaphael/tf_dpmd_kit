import glob
import os

list2d_snap = [
    #[100,150],
    #range(140,157)
]

for list_snap in list2d_snap:
    for int_snap in list_snap:
        str_dir = f'snap_{int_snap:0>5d}'
        if os.path.exists(str_dir):
            print(str_dir)
            os.rename(str_dir, 'bk.'+str_dir)
        str_file = f'{str_dir}.POSCAR'
        if os.path.isfile(str_file):
            print(str_file)
            os.rename(str_file, 'bk.'+str_file)
