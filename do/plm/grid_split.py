import os

str_save = 'fes.dist_vp_c.grid'

if os.path.isfile(str_save):
   os.rename(str_save, f'bk.{str_save}')

with open('../fes.dist_vp_c.grid', 'r') as fp_r:
    int_i = 0
    str_write = fp_r.readline()
    for str_line in fp_r:

        str_write += str_line
        for str_line in fp_r:
            if str_line[0] == '#':
                str_write += str_line
            else:
                break

        str_write += str_line
        for str_line in fp_r:
            if str_line[0] != '#':
                str_write += str_line
            else:
                break

        if os.path.isfile(str_save):
            os.rename(str_save, f'analysis.{int_i}.{str_save}')
            int_i += 1
        with open(str_save, 'w') as fp_w:
            fp_w.write(str_write)

        str_write = str_line
