import ase.io

dict_pwscfin = {}
dict_pwscfin['CONTROL'] = {
    'dt': 20.0,
    'nstep': 80000,
    'max_seconds': 82800,
    'disk_io': 'none',
    }
dict_pwscfin['SYSTEM'] = {
    'input_dft': 'scan',
    'ecutwfc': 150,
    }
dict_pwscfin['ELECTRONS'] = {
    'electron_dynamics': 'cp-bo',
    'emass': 1000.0,
    'emass_emin': 200,
    'emass_cutoff': 6.0,
    'electron_damping_emin': 0.7,
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
    'O': 'O_ONCV_PBE-1.2.upf',
    'H': 'H_ONCV_PBE-1.2.upf'
    }

atoms_poscar = ase.io.read(
    filename = 'POSCAR',
    format = 'vasp',
    )

ase.io.write(
    filename = 'pwscf.in',
    images = atoms_poscar,
    format = 'espresso-in',
    input_data = dict_pwscfin,
    pseudopotentials = dict_pwpseudop )

