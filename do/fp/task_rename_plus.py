import os
import glob

def run():

    list_task = glob.glob('task.*.*')
    for str_task in list_task:
        list_tmp = str_task.split('.')
        str_new = list_tmp[0]+'.'+list_tmp[1]+'.1'+list_tmp[2]
        print(str_task, 'to', str_new)
        os.rename(str_task, str_new)
 
run()
