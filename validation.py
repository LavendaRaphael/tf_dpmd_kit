import dpdata
import numpy as np
from string import Template
from dpgen.generator.lib.utils import create_path
import os
import ase.io
from dpgen.remote.decide_machine import convert_mdata
import json
from dpgen.dispatcher.Dispatcher import make_submission
import glob

fp_name = '02.fp'

def def_id2dir(int_id):
    return f"task.{int_id:06d}"

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
def make_fp (
             jdata,
             mdata) :
    fp_style = jdata['fp_style']

    if fp_style == "pwscf" :
        make_fp_pwscf(jdata)
    else :
        raise RuntimeError ("unsupported fp style")
    # Copy user defined forward_files
    work_path = fp_name
    #symlink_user_forward_files(mdata=mdata, task_type="fp", work_path=work_path)

def run_fp_inner (
                  jdata,
                  mdata,
                  forward_files,
                  backward_files,
                  log_file = "fp.log",
                  forward_common_files=[]) :
    fp_command = mdata['fp_command']
    fp_group_size = mdata['fp_group_size']
    fp_resources = mdata['fp_resources']

    work_path = fp_name

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return

    fp_run_tasks = fp_tasks
    run_tasks = [os.path.basename(ii) for ii in fp_run_tasks]

    user_forward_files = mdata.get("fp" + "_user_forward_files", [])
    forward_files += [os.path.basename(file) for file in user_forward_files]
    backward_files += mdata.get("fp" + "_user_backward_files", [])

    submission = make_submission(
        mdata['fp_machine'],
        mdata['fp_resources'],
        commands=[fp_command],
        work_path=work_path,
        run_tasks=run_tasks,
        group_size=fp_group_size,
        forward_common_files=forward_common_files,
        forward_files=forward_files,
        backward_files=backward_files,
        outlog = log_file,
        errlog = log_file)
    submission.run_submission()

def run_fp (
            jdata,
            mdata) :
    '''
    https://github.com/deepmodeling/dpgen/blob/dc57c9a2e7ada3cf1c48c6b8bc9d475bb95ebf40/dpgen/generator/run.py#L2426-L2470
    '''
    fp_style = jdata['fp_style']
    fp_pp_files = jdata['fp_pp_files']
    if fp_style == "pwscf" :
        forward_files = ['input'] #+ fp_pp_files
        backward_files = ['output']
        run_fp_inner(jdata, mdata,  forward_files, backward_files, log_file = 'output')
    else :
        raise RuntimeError ("unsupported fp style")

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
