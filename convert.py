import ase.io

def def_poscar2vccpbo():
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
        'nhpcl': 4,
        'nhptyp': 2,
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


def def_poscar2cpgs():

    dict_pwscfin = {}
    dict_pwscfin['CONTROL'] = {
        'restart_mode': 'from_scratch',
        'dt': 2.0,
        'nstep': 20000,
        'max_seconds': 82800,
        'tstress': True,
        'tprnfor': True
        }
    dict_pwscfin['SYSTEM'] = {
        'input_dft': 'scan',
        'ecutwfc': 150,
        'ibrav': 1
        }
    dict_pwscfin['ELECTRONS'] = {
        'electron_dynamics': 'damp',
        'electron_damping': 0.2,
        'emass': 100.0,
        'emass_cutoff': 25.0,
        'ortho_max': 10000,
        'electron_maxstep': 500
        }
    dict_pwscfin['IONS'] = {
        'ion_dynamics': 'none',
        'ion_radius(1)': 1.4,
        'ion_radius(2)': 1.4,
        'ion_radius(3)': 1.4
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

