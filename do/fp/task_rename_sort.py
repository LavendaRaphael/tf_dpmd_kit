import os
import glob

def run():

    list_task = glob.glob('task.*')
    list_system = []
    for str_task in list_task :
        list_system.append(str_task.split('.')[1])
    list_system = list(set(list_system))

    for str_system in list_system:
        list_task = glob.glob(f'task.{str_system}.*')
        int_count = 0
        for int_i in range(len(list_task)):
            str_task = f'task.{str_system}.{int_i:0>6d}'
            if not os.path.exists(str_task):
                int_count += 1
                str_rename = list_task[-int_count]
                print(str_rename, 'to', str_task)
                os.rename(str_rename, str_task)
        
run()
