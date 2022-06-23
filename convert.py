import dpdata
import ase.io
import numpy
import os
import json

def def_cp2pwscf(
        np_snap
        ):
    print(np_snap)
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
    
    str_in = 'pwscf.in'
    str_log = 'pwscf.out'
    str_command = "if [ -f ../please.stop ]; then true; else"
    str_command += " cat $PBS_NODEFILE|sort -u|xargs echo 'NODE:'"
    str_command += " && mpirun qe.7.0_libxc_pw.x < "+str_in
    str_command += " && if [ -f ../please.continue ]; then false; fi; fi"
    dict_task = {
        "command": str_command,
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

    if (not os.path.exists('snap')):
        os.mkdir('snap')
    os.chdir('snap')
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

def def_dump2ase(
        array_id,
        type_map,       # ["O", "H"]
        ):
    dp_sys = dpdata.System()
    for int_i in array_id:
        str_origin = f"traj/{int_i}.lammpstrj"
        dp_tmp = dpdata.System(str_origin, fmt='lammps/dump', type_map=type_map)

        dp_sys.append(dp_tmp)
    print(dp_sys)
    return dp_sys.to('ase/structure')

def def_poscar2cp_cp():
    dict_pwscfin = {}
    dict_pwscfin['CONTROL'] = {
        'calculation': 'cp',
        'dt': 2.0,
        'nstep': 20000,
        'max_seconds': 82800,
        'isave': 500,
        'iprint': 10,
        'tstress': True,
        'tprnfor': True,
        }
    dict_pwscfin['SYSTEM'] = {
        'input_dft': 'scan',
        'ecutwfc': 85.0,
        'ibrav': 1
        }
    dict_pwscfin['ELECTRONS'] = {
        'electron_dynamics': 'verlet',
        'emass': 100.0,
        'emass_cutoff': 25.0,
        'ortho_max': 800
        }
    dict_pwscfin['IONS'] = {
        'ion_dynamics': 'verlet',
        'ion_temperature': 'nose',
        'tempw': 330,
        'fnosep': 60,
        'nhptyp': 0,
        'ndega': -3,
        }

    dict_pwpseudop = {
        'O': 'O_HSCV_PBE-1.0.UPF',
        'H': 'H_HSCV_PBE-1.0.UPF',
        'C': 'C_HSCV_PBE-1.0.UPF'
        }

    atoms_poscar = ase.io.read(
        filename = 'POSCAR',
        format = 'vasp',
        )

    atoms_poscar.set_masses( [ 2.0141 if x==1.008 else x for x in atoms_poscar.get_masses() ] )

    ase.io.write(
        filename = 'cp.in',
        images = atoms_poscar,
        format = 'espresso-in',
        input_data = dict_pwscfin,
        pseudopotentials = dict_pwpseudop )

def def_poscar2cp_cpbo():
    dict_pwscfin = {}
    dict_pwscfin['CONTROL'] = {
        'calculation': 'cp',
        'dt': 20.0,
        'nstep': 20000,
        'max_seconds': 82800,
        'isave': 50,
        'iprint': 5,
        'tstress': True,
        'tprnfor': True,
        }
    dict_pwscfin['SYSTEM'] = {
        'input_dft': 'scan',
        'ecutwfc': 85.0,
        'ibrav': 1
        }
    dict_pwscfin['ELECTRONS'] = {
        'electron_dynamics': 'cp-bo',
        'emass_emin': 200,
        'emass_cutoff_emin': 6.0,
        'electron_damping_emin': 0.7,
        'dt_emin': 5.0,
        'emass': 1000.0,
        'emass_cutoff': 6.0,
        'ortho_max': 300,
        'electron_maxstep': 500
        }
    dict_pwscfin['IONS'] = {
        'ion_dynamics': 'verlet',
        'ion_temperature': 'nose',
        'fnosep': 60,
        'tempw': 330,
        'nhpcl': 4,
        'nhptyp': 2,
        'ndega': -3,
        }

    dict_pwpseudop = {
        'O': 'O_HSCV_PBE-1.0.UPF',
        'H': 'H_HSCV_PBE-1.0.UPF',
        'C': 'C_HSCV_PBE-1.0.UPF'
        }

    atoms_poscar = ase.io.read(
        filename = 'POSCAR',
        format = 'vasp',
        )

    atoms_poscar.set_masses( [ 2.0141 if x==1.008 else x for x in atoms_poscar.get_masses() ] )

    ase.io.write(
        filename = 'cp.in',
        images = atoms_poscar,
        format = 'espresso-in',
        input_data = dict_pwscfin,
        pseudopotentials = dict_pwpseudop )

def def_poscar2cp_gs():

    dict_pwscfin = {}
    dict_pwscfin['CONTROL'] = {
        'calculation': 'cp',
        'restart_mode': 'from_scratch',
        'dt': 2.0,
        'nstep': 20000,
        'max_seconds': 82800,
        'tstress': True,
        'tprnfor': True
        }
    dict_pwscfin['SYSTEM'] = {
        'input_dft': 'scan',
        'ecutwfc': 85.0,
        'ibrav': 1,
        }
    dict_pwscfin['ELECTRONS'] = {
        'electron_dynamics': 'damp',
        'electron_damping': 0.20,
        'emass': 100.0,
        'emass_cutoff': 25.0,
        'ortho_max': 800
        }
    dict_pwscfin['IONS'] = {
        'ion_dynamics': 'none',
        'ion_radius(1)': 1.4,
        'ion_radius(2)': 1.4,
        'ion_radius(3)': 1.4,
        }
    
    dict_pwpseudop = {
        'O': 'O_HSCV_PBE-1.0.UPF',
        'H': 'H_HSCV_PBE-1.0.UPF',
        'C': 'C_HSCV_PBE-1.0.UPF'
        }
    
    atoms_poscar = ase.io.read(
        filename = 'POSCAR',
        format = 'vasp',
        )
    
    atoms_poscar.set_masses( [ 2.0141 if x==1.008 else x for x in atoms_poscar.get_masses() ] )
    
    ase.io.write(
        filename = 'gs.in',
        images = atoms_poscar,
        format = 'espresso-in',
        input_data = dict_pwscfin,
        pseudopotentials = dict_pwpseudop )

def def_poscar2vccp_cpbo():
    dict_pwscfin = {}
    dict_pwscfin['CONTROL'] = {
        'calculation': 'vc-cp',
        'dt': 20.0,
        'nstep': 20000,
        'max_seconds': 82800,
        'isave': 50,
        'iprint': 5,
        'tstress': True,
        'tprnfor': True,
        }
    dict_pwscfin['SYSTEM'] = {
        'input_dft': 'scan',
        'ecutwfc': 150,
        'ecfixed': 130,
        'qcutz': 200,
        'q2sigma': 15,
        'ibrav': 1
        }
    dict_pwscfin['ELECTRONS'] = {
        'electron_dynamics': 'cp-bo',
        'emass_emin': 200,
        'emass_cutoff_emin': 6.0,
        'electron_damping_emin': 0.7,
        'dt_emin': 5.0,
        'emass': 1000.0,
        'emass_cutoff': 6.0,
        'ortho_max': 300,
        'electron_maxstep': 500
        }
    dict_pwscfin['IONS'] = {
        'ion_dynamics': 'verlet',
        'ion_temperature': 'nose',
        'fnosep': 60,
        'tempw': 330,
        'nhptyp': 0,
        'ndega': -3,
        }
    dict_pwscfin['CELL'] = {
        'cell_dynamics': 'pr',
        'press': 0.001,
        'cell_dofree': 'volume',
        }

    dict_pwpseudop = {
        'O': 'O_HSCV_PBE-1.0.UPF',
        'H': 'H_HSCV_PBE-1.0.UPF',
        'C': 'C_HSCV_PBE-1.0.UPF'
        }

    atoms_poscar = ase.io.read(
        filename = 'POSCAR',
        format = 'vasp',
        )

    atoms_poscar.set_masses( [ 2.0141 if x==1.008 else x for x in atoms_poscar.get_masses() ] )

    ase.io.write(
        filename = 'cp.in',
        images = atoms_poscar,
        format = 'espresso-in',
        input_data = dict_pwscfin,
        pseudopotentials = dict_pwpseudop )

    float_cell = atoms_poscar.get_cell()[0,0]
    float_cell_ref = float_cell * 1.03
    with open('cp.in', 'a') as open_in:
        open_in.write(f'REF_CELL_PARAMETERS angstrom\n  {float_cell_ref} 0.0 0.0\n  0.0 {float_cell_ref} 0.0\n  0.0 0.0 {float_cell_ref}')

def def_poscar2vccp_cp():
    dict_pwscfin = {}
    dict_pwscfin['CONTROL'] = {
        'calculation': 'vc-cp',
        'dt': 2.0,
        'nstep': 20000,
        'max_seconds': 82800,
        'isave': 500,
        'iprint': 10,
        'tstress': True,
        'tprnfor': True,
        }
    dict_pwscfin['SYSTEM'] = {
        'input_dft': 'scan',
        'ecutwfc': 150,
        'ecfixed': 130,
        'qcutz': 200,
        'q2sigma': 15,
        'ibrav': 1
        }
    dict_pwscfin['ELECTRONS'] = {
        'electron_dynamics': 'verlet',
        'emass': 100.0,
        'emass_cutoff': 25.0,
        'ortho_max': 800
        }
    dict_pwscfin['IONS'] = {
        'ion_dynamics': 'verlet',
        'ion_temperature': 'nose',
        'fnosep': 60,
        'tempw': 330,
        'nhptyp': 0,
        'ndega': -3,
        }
    dict_pwscfin['CELL'] = {
        'cell_dynamics': 'pr',
        'press': 0.001,
        'cell_dofree': 'volume',
        }

    dict_pwpseudop = {
        'O': 'O_HSCV_PBE-1.0.UPF',
        'H': 'H_HSCV_PBE-1.0.UPF',
        'C': 'C_HSCV_PBE-1.0.UPF'
        }

    atoms_poscar = ase.io.read(
        filename = 'POSCAR',
        format = 'vasp',
        )

    atoms_poscar.set_masses( [ 2.0141 if x==1.008 else x for x in atoms_poscar.get_masses() ] )

    ase.io.write(
        filename = 'cp.in',
        images = atoms_poscar,
        format = 'espresso-in',
        input_data = dict_pwscfin,
        pseudopotentials = dict_pwpseudop )

    float_cell = atoms_poscar.get_cell()[0,0]
    float_cell_ref = float_cell * 1.03
    with open('cp.in', 'a') as open_in:
        open_in.write(f'REF_CELL_PARAMETERS angstrom\n  {float_cell_ref} 0.0 0.0\n  0.0 {float_cell_ref} 0.0\n  0.0 0.0 {float_cell_ref}')

def def_poscar2vccp_gs():

    dict_pwscfin = {}
    dict_pwscfin['CONTROL'] = {
        'calculation': 'vc-cp',
        'restart_mode': 'from_scratch',
        'dt': 2.0,
        'nstep': 20000,
        'max_seconds': 82800,
        'tstress': True,
        'tprnfor': True
        }
    dict_pwscfin['SYSTEM'] = {
        'input_dft': 'scan',
        'ecutwfc': 150.0,
        'ibrav': 1,
        'ecfixed': 130.0,
        'qcutz': 200.0,
        'q2sigma': 15.0,
        }
    dict_pwscfin['ELECTRONS'] = {
        'electron_dynamics': 'damp',
        'emass': 100.0,
        'emass_cutoff': 25.0,
        'ortho_max': 1000,
        }
    dict_pwscfin['IONS'] = {
        'ion_dynamics': 'none',
        }
    dict_pwscfin['CELL'] = {
        'cell_dynamics': 'none'
        }

    dict_pwpseudop = {
        'O': 'O_HSCV_PBE-1.0.UPF',
        'H': 'H_HSCV_PBE-1.0.UPF',
        'C': 'C_HSCV_PBE-1.0.UPF'
        }
    
    atoms_poscar = ase.io.read(
        filename = 'POSCAR',
        format = 'vasp',
        )
    
    atoms_poscar.set_masses( [ 2.0141 if x==1.008 else x for x in atoms_poscar.get_masses() ] )
    
    ase.io.write(
        filename = 'gs.in',
        images = atoms_poscar,
        format = 'espresso-in',
        input_data = dict_pwscfin,
        pseudopotentials = dict_pwpseudop )
    
    float_cell = atoms_poscar.get_cell()[0,0]
    float_cell_ref = float_cell * 1.03
    with open('gs.in', 'a') as open_in:
        open_in.write(f'REF_CELL_PARAMETERS angstrom\n  {float_cell_ref} 0.0 0.0\n  0.0 {float_cell_ref} 0.0\n  0.0 0.0 {float_cell_ref}')
