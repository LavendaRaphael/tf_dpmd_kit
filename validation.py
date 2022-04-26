import dpdata
import numpy as np
from dpgen.generator.lib.utils import create_path
import os
import ase.io
from dpgen.remote.decide_machine import convert_mdata
import json

def make_fp_pwscf():
    array_id = np.arange(
        start = 0,
        stop = 5002,
        step = 20,
        dtype = int,
        )
    
    dp_sys = dpdata.System(
        file_name = '01.dpmd_equil/traj.dump',
        fmt = "lammps/dump",
        type_map = ["O", "H"]
        )
    print(dp_sys)
    atoms = dp_sys.to('ase/structure')
    
    dict_pwscfin = {}
    dict_pwscfin['CONTROL'] = {
        'tstress': True,
        'tprnfor': True,
        'disk_io': 'none',
        }
    dict_pwscfin['SYSTEM'] = {
        'input_dft': 'scan',
        'ecutwfc': 150,
        'nosym': True,
        }
    dict_pwscfin['ELECTRONS'] = {
        'electron_maxstep': 500
        }
    dict_pwpseudop = {
        'O': 'O_HSCV_PBE-1.0.UPF',
        'H': 'H_HSCV_PBE-1.0.UPF'
        }
    
    create_path( fp_name )
    os.chdir(fp_name)
    
    str_dir_cwd = os.getcwd()
    for int_id in array_id:
        print(int_id)
        str_subdir = def_id2dir(int_id)
    
        create_path( str_subdir )
        os.chdir(str_subdir)
    
        ase.io.write(
            filename = 'input',
            images = atoms[int_id],
            format = 'espresso-in',
            input_data = dict_pwscfin,
            pseudopotentials = dict_pwpseudop )
        
        os.chdir(str_dir_cwd)

def read_para_machine(
        param_file = 'param.json', 
        machine_file = 'machine.json'
        ) :
    '''
    https://github.com/deepmodeling/dpgen/blob/dc57c9a2e7ada3cf1c48c6b8bc9d475bb95ebf40/dpgen/generator/run.py#L2911-L2939
    '''
    
    with open (param_file, 'r') as fp :
        jdata = json.load (fp)
    with open (machine_file, 'r') as fp:
        mdata = json.load (fp)

    mdata = convert_mdata(mdata)

    return jdata, mdata
