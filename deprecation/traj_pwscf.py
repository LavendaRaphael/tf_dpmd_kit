#!/usr/bin/env python
# coding: utf-8

import ase.io
import os
import shutil
import numpy
import group_module
import dpdata
import glob
import os

def def_pwscf2dpnpy_merge(
        class_paras,
        int_copy,
        str_outdir,
        ):
    
    dp_sys = dpdata.MultiSystems.from_dir(
        dir_name = class_paras.str_workdir,
        file_name = 'pwscf.out',
        fmt = 'qe/pw/scf'
        )[0]
    print(dp_sys)

    dp_syscopy = dp_sys.copy()
    for int_i in range(int_copy-1):
        dp_syscopy.append(dp_sys)
    print(dp_syscopy)

    dp_syscopy.to(
        'deepmd/npy',
        str_outdir,
        )

    os.chdir( str_cwd )

def def_allrename(
        dir_name,
        str_origin,
        str_new,
        ):

    list_file = sorted(glob.glob('./{}/**/{}'.format(dir_name, str_origin), recursive=True))
    for str_file in list_file:
        print(str_file)
        os.rename(str_file, str_file.replace(str_origin, str_new))

def def_inputinit_pwscf(
        class_paras,
        ):
    atoms_traj = []
    
    str_cwd = os.getcwd()
    os.chdir( class_paras.str_workdir )

    for int_id in class_paras.array_id:
        print(int_id)
        str_subdir = class_paras.def_id2dir( int_id )

        os.chdir(str_subdir)
        atoms_poscar = ase.io.read(
            filename = 'POSCAR',
            format = 'vasp',
            )

        if (not os.path.isdir(  class_paras.str_compsubdir )):
            os.mkdir( class_paras.str_compsubdir )
        os.chdir( class_paras.str_compsubdir )
        ase.io.write(
            filename = 'pwscf.in',
            images = atoms_poscar, 
            format = 'espresso-in', 
            input_data = class_paras.dict_pwscfin, 
            pseudopotentials = class_paras.dict_pwpseudop )

        os.chdir('..')
        os.chdir('..')
    os.chdir( str_cwd )

def def_traj2poscar(
        class_paras,
        ):
    atoms_traj = []
    
    for list_file in class_paras.list2d_files:
        str_filename = list_file[0]
        str_format = list_file[1]
        str_ppcode = list_file[2]
        if (str_ppcode == 'ase'):
            atoms_tmp = ase.io.read(
                filename = str_filename, 
                format = str_format, 
                index = ':')
        elif (str_ppcode == 'dpdata'):
            atoms_tmp = def_dpsys2ase(
                str_filename = str_filename,
                str_format = str_format,
                )
        atoms_traj.extend( atoms_tmp )
    int_nframes = len(atoms_traj)
    print('int_nframes = ', int_nframes)
   
    str_cwd = os.getcwd()
    if (not os.path.isdir(  class_paras.str_workdir )):
        os.mkdir( class_paras.str_workdir )
    os.chdir( class_paras.str_workdir )

    for int_id in class_paras.array_id:
        print(int_id)
        str_subdir = class_paras.def_id2dir( int_id )
        if (not os.path.isdir(str_subdir)):
            #shutil.rmtree(str_subdir)
            os.mkdir(str_subdir)
        os.chdir(str_subdir)

        ase.io.write(
            filename = 'POSCAR',
            images = atoms_traj[ int_id ], 
            format = 'vasp', 
            )

        os.chdir('..')
    os.chdir( str_cwd )

def def_dpsys2ase( str_filename, str_format ):
    
    dp_sys = dpdata.System(
        file_name = str_filename,
        fmt = str_format,
        )
    return dp_sys.to('ase/structure')

class class_paras( group_module.class_subparas ):
    def __init__(self):
        super( class_paras, self ).__init__()
    
        self._dict_pwscfin = {}
        self._dict_pwscfin['CONTROL'] = {
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none',
            }
        self._dict_pwscfin['SYSTEM'] = {
            'input_dft': 'scan',
            'ecutwfc': 150,
            'nosym': True,
            }
        self._dict_pwscfin['ELECTRONS'] = {
            'electron_maxstep': 500
            }
        #
        self._dict_pwpseudop = {
            'O': 'O_HSCV_PBE-1.0.UPF',
            'H': 'H_HSCV_PBE-1.0.UPF'
            }

    @property
    def dict_pwscfin(self):
        return self._dict_pwscfin

    @property
    def dict_pwpseudop(self):
        return self._dict_pwpseudop

    @property
    def list2d_files(self):
        return self._list2d_files
    @list2d_files.setter
    def list2d_files(self, value):
        self._list2d_files = value
    
