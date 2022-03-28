#!/usr/bin/env python
# coding: utf-8

import ase.io
import os
import shutil
import numpy
import group_module

def def_inputinit_pwscf(
        class_paras,
        ):
    atoms_traj = []
    for str_filename in class_paras.list1d_filename:
        atoms_tmp = ase.io.read(
            filename = str_filename, 
            format = class_paras.str_format, 
            index = ':')
        atoms_traj.extend( atoms_tmp )
    int_nframes = len(atoms_traj)
    print('int_nframes = ', int_nframes)
    
    if (not os.path.isdir(  class_paras.str_workdir )):
        os.mkdir( class_paras.str_workdir )
    os.chdir( class_paras.str_workdir )

    for int_id in class_paras.array_id:
        print(int_id)
        str_subdir = class_paras.def_id2dir( int_id )
        if os.path.isdir(str_subdir):
            shutil.rmtree(str_subdir)
        os.mkdir(str_subdir)
        os.chdir(str_subdir)

        ase.io.write(
            filename = 'pwscf.in',
            images = atoms_traj[ int_id ], 
            format = 'espresso-in', 
            input_data = class_paras.dict_pwscfin, 
            pseudopotentials = class_paras.dict_pwpseudop )

        os.chdir('..')

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
    def list1d_filename(self):
        return self._list1d_filename
    @list1d_filename.setter
    def list1d_filename(self, value):
        self._list1d_filename = value
    
    @property
    def str_format(self):
        return self._str_format
    @str_format.setter
    def str_format(self, value):
        self._str_format = value

