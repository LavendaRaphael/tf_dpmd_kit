import glob
import os

dir_name = './'
str_origin = 'pwscf.log'
str_new = 'pwscf.out'

list_file = sorted(glob.glob('./{}/**/{}'.format(dir_name, str_origin), recursive=True))
for str_file in list_file:
    print(str_file)
    os.rename(str_file, str_file.replace(str_origin, str_new))

