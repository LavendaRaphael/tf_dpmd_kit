#!/usr/bin/env python
# coding: utf-8

import ase.io
import os

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
    print(int_nframes)
    
    array_select = numpy.linspace( 
        0, 
        int_nframes, 
        num = class_paras.int_select, 
        dtype = int, 
        endpoint = False )

    os.chdir( class_paras.str_dir )
    for int_i in array_select:
        str_i = "snap_%05d" % int_i
        print(str_i)
        os.mkdir(str_i)

        ase.io.write(
            filename = 'pwscf.in',
            images = atoms_i, 
            format = 'espresso-in', 
            input_data = class_paras.dict_pwscfin, 
            pseudopotentials = class_paras.dict_pwpseudop )

class class_paras():
    def __init__(self):

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
    def str_dir(self):
        return self._str_dir
    @str_dir.setter
    def str_dir(self, value):
        self._str_dir = value

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

    @property
    def int_select(self):
        return self._int_select
    @int_select.setter
    def int_select(self, value):
        self._int_select = value


