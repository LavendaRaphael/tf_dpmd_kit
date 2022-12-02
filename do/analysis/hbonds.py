import MDAnalysis as mda
import numpy as np
from tf_dpmd_kit import analysis
import os

# read structrue

#TT
dict_sel = {
    'o_w': 'index 0:125',
    'o_0_2': 'index 126 128',
    'o_1': 'index 127',
    'h_w': 'index 129:380',
    'h_0_1': 'index 381 382',
    'c': 'index 383'
}

dict_atompair = {
    #'h_w.o_1': ('h_w', 'o_1'),
    #'h_w.o_0_2': ('h_w', 'o_0_2'),
    #'h_0_1.o_w': ('h_0_1','o_w'),
    'h_w.o_w': ('h_w','o_w')
}

dict_ave = {
    'h_w.o_w': 126.0/2.0
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
    dict_atompair: list,
    dict_donors: dict = {},
    dict_ave: dict = {}
)->None:

    mda_u = mda.Universe(str_topology, list_traj, topology_format="DATA", format="LAMMPSDUMP")
    mda_u.select_atoms("type 1").types = 'O'
    mda_u.select_atoms("type 2").types = 'H'
    mda_u.select_atoms("type 3").types = 'C'
    print(mda_u.trajectory)

    for str_label, str_sel in dict_sel.items():
        print(str_label, mda_u.select_atoms(str_sel).indices)

    for str_label, tup_atompair in dict_atompair.items():

        if str_label in dict_donors:
            donors_sel = dict_donors[str_label]
            update_selections = True
        else:
            donors_sel = None
            update_selections = False

        if str_label in dict_ave:
            float_ave = dict_ave[str_label]
        else:
            float_ave = 1.0

        analysis.hbonds(
            universe = mda_u,
            hydrogens_sel = dict_sel[tup_atompair[0]],
            acceptors_sel = dict_sel[tup_atompair[1]],
            donors_sel = donors_sel,
            update_selections = update_selections,
            str_save = f'hbonds.{tup_atompair[0]}.{tup_atompair[1]}.csv',
            d_a_cutoff = 3.5,
            d_h_a_angle_cutoff = 150,
            float_ave = float_ave,
        )

str_cwd = os.getcwd()
for str_dir in list_dir:
    print(str_dir)
    os.chdir(str_dir)
    run(
        str_topology = '../../lmp.data',
        list_traj = ['traj.0.lammpstrj', 'traj.lammpstrj'],
        dict_sel = dict_sel,
        dict_atompair = dict_atompair,
        dict_ave = dict_ave,
    )
    os.chdir(str_cwd)
