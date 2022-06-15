import dpdata
import ase.io
import numpy
import os
import json

np_snap = numpy.arange(
    start = 0,
    stop = 2000,
    step = 80,
    )

dp_sys = dpdata.System(
    file_name = 'cp',
    fmt = 'qe/cp/traj',
    )
print(dp_sys)
print(dp_sys['cells'][0])

dict_pwscfin = {}
dict_pwscfin['CONTROL'] = {
    'tstress': True,
    'tprnfor': True,
    'disk_io': 'none',
    }
dict_pwscfin['SYSTEM'] = {
    'input_dft': 'scan',
    'ecutwfc': 150,
    }

dict_pwpseudop = {
    'O': 'O_HSCV_PBE-1.0.UPF',
    'H': 'H_HSCV_PBE-1.0.UPF',
    'C': 'C_HSCV_PBE-1.0.UPF'
    }

if (not os.path.exists('snap')):
    os.mkdir('snap')
os.chdir('snap')

str_in = 'pwscf.in'
str_log = 'pwscf.log'
dict_task = {
    "command": "mpirun qe.7.0_libxc_pw.x < "+str_in,
    "forward_files": [
        str_in
        ],
    "backward_files": [
        str_log
        ],
    'outlog': str_log,
    'errlog': str_log
    }

ase_atoms = dp_sys.to('ase/structure')
for int_snap in np_snap:
    str_dir = f'snap_{int_snap:0>5d}'
    if (not os.path.exists(str_dir)):
        os.mkdir(str_dir)
    os.chdir(str_dir)

    ase.io.write(
        filename = str_in,
        images = ase_atoms[int_snap],
        format = 'espresso-in',
        input_data = dict_pwscfin,
        pseudopotentials = dict_pwpseudop 
        )

    dict_task['task_work_path'] = str_dir
    with open('task.json', 'w') as open_json:
        json.dump( dict_task, fp=open_json, indent=4 )

    os.chdir('..')
os.chdir('..')


