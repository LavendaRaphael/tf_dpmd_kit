import MDAnalysis as mda
import numpy as np
from tf_dpmd_kit import analysis
import os

# read structrue

dict_sel = {
    'o_w': 'index 0:125',
    'o_0': 'index 126',
    'o_1_2': 'index 127 128',
    'h_w': 'index 129:380 382',
    'h_0': 'index 381',
    'c': 'index 383',
}
dict_atompair = {
    'h_w.o_0': ('h_w', 'o_0'),
    'h_w.o_1_2': ('h_w', 'o_1_2'),
    'h_0.o_w': ('h_0', 'o_w'),
}
dict_donors = {
    'h_w': 'o_w',
}

list_dir = [
    '280K/hbonds',
    '290K/hbonds',
    '300K/hbonds',
    '310K/hbonds',
    '320K/hbonds',
]

def run(
    str_topology: str,
    list_traj: list,
    dict_sel: dict,
    dict_atompair: dict,
    dict_donors: dict,
)->None:

    mda_u = mda.Universe(str_topology, list_traj, topology_format="DATA", format="LAMMPSDUMP")
    mda_u.select_atoms("type 1").types = 'O'
    mda_u.select_atoms("type 2").types = 'H'
    mda_u.select_atoms("type 3").types = 'C'
    print(mda_u.trajectory)

    for str_label, str_sel in dict_sel.items():
        print(str_label, mda_u.select_atoms(str_sel).indices)

    for str_label, tup_atompair in dict_atompair.items():

        hydrogens_sel = dict_sel[tup_atompair[0]]
        acceptors_sel = dict_sel[tup_atompair[1]]

        if hydrogens_sel in dict_donors:
            donors_sel = dict_sel[dict_donors[hydrogens_sel]]
            update_selections = True
        else:
            donors_sel = None
            update_selections = False

        analysis.hbonds(
            universe = mda_u,
            hydrogens_sel = hydrogens_sel,
            acceptors_sel = acceptors_sel,
            donors_sel = donors_sel,
            update_selections = update_selections,
            str_save = f'hbonds.{tup_atompair[0]}.{tup_atompair[1]}.csv',
            d_a_cutoff = 3.5,
            d_h_a_angle_cutoff = 150,
        )

str_cwd = os.getcwd()
for str_dir in list_dir:
    print(str_dir)
    os.chdir(str_dir)
    run(
        str_topology = '../../lmp.data',
        list_traj = ['traj.lammpstrj.0', 'traj.lammpstrj'],
        dict_sel = dict_sel,
        dict_atompair = dict_atompair,
        dict_donors = dict_donors,
    )
    os.chdir(str_cwd)
