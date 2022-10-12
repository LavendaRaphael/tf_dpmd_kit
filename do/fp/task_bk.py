import glob
import os

list2d_task = [
    [71, [4]],
]

for list_task in list2d_task:
    int_task = list_task[0]
    for int_snap in list_task[1]:
        str_dir = f'task.{int_task:0>3d}.{int_snap:0>6d}'
        if os.path.exists(str_dir):
            print(str_dir)
            os.rename(str_dir, 'bk.'+str_dir)
        str_file = f'{str_dir}.POSCAR'
        if os.path.isfile(str_file):
            print(str_file)
            os.rename(str_file, 'bk.'+str_file)
