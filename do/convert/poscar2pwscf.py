import ase.io

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

