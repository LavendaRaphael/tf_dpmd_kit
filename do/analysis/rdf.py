import MDAnalysis as mda
import numpy as np
from tf_dpmd_kit import analysis
import os

#CC
dict_sel = {
    'o_w': 'index 0:125',
    'o_0_2': 'index 126 128',
    'o_1': 'index 127',
    'h_w': 'index 129:380',
    'h_0_1': 'index 381 382',
    'c': 'index 383'
}
list_atompair = [
    ('o_1', 'h_w'),
    ('o_0_2','h_w'),
    ('h_0_1','o_w'),
    ('o_w','o_w')
]

list_snaprange = [
    (     0, 100000),
    (100000, 200000),
    (200000, 300000),
    (300000, 400000),
    (100000, 400000),
]

list_dir = [
    '280K/rdf',
    '290K/rdf',
    '300K/rdf',
    '310K/rdf',
    '320K/rdf',
]

def run(
    str_topology: str,
    list_traj: list,
    dict_sel: dict,
    list_snaprange: list,
    list_atompair: list,
)->None:

    mda_u = mda.Universe(str_topology, list_traj, topology_format="LAMMPSDUMP", format="LAMMPSDUMP")
    mda_u.select_atoms("type 1").types = 'O'
    mda_u.select_atoms("type 2").types = 'H'
    mda_u.select_atoms("type 3").types = 'C'
    print(mda_u.trajectory)
    
    for str_label, str_sel in dict_sel.items():
        print(str_label, mda_u.select_atoms(str_sel).indices)
    
    for tup_snaprange in list_snaprange:
        for tup_atompair in list_atompair:
            analysis.rdf(
                mda_atomgroup_0 = mda_u.select_atoms(dict_sel[tup_atompair[0]]),
                mda_atomgroup_1 = mda_u.select_atoms(dict_sel[tup_atompair[1]]),
                tup_snaprange = tup_snaprange,
                str_save = f'rdf.{tup_atompair[0]}.{tup_atompair[1]}.{tup_snaprange[0]:07d}_{tup_snaprange[1]:07d}.csv',
                tup_rrange = (1.0,6.0)
            )

str_cwd = os.getcwd()
for str_dir in list_dir:
    print(str_dir)
    os.chdir(str_dir)
    run(
        str_topology = 'traj.lammpstrj',
        list_traj = ['traj.0.lammpstrj', 'traj.lammpstrj'],
        dict_sel = dict_sel,
        list_snaprange = list_snaprange,
        list_atompair = list_atompair,
    )
    os.chdir(str_cwd)
