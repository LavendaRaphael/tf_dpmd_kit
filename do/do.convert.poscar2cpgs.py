import ase.io

dict_pwscfin = {}
dict_pwscfin['CONTROL'] = {
    'dt': 2.0,
    'nstep': 20000,
    'max_seconds': 82800,
    }
dict_pwscfin['SYSTEM'] = {
    'input_dft': 'scan',
    'ecutwfc': 150,
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

