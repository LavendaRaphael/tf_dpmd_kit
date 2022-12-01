import MDAnalysis as mda
import numpy as np
from tf_dpmd_kit import analysis

# read structrue

mda_u = mda.Universe('lmp.data', ['traj.lammpstrj'], topology_format="DATA", format="LAMMPSDUMP")
mda_u.select_atoms("type 1").types = 'O'
mda_u.select_atoms("type 2").types = 'H'
mda_u.select_atoms("type 3").types = 'C'
print(mda_u.trajectory)

#TT
dict_sel = {
    'o_w': 'index 0:125',
    'o_0_2': 'index 126 128',
    'o_1': 'index 127',
    'h_w': 'index 129:380',
    'h_0_1': 'index 381 382',
    'c': 'index 383'
}

for str_label, str_sel in dict_sel.items():
    print(str_label, mda_u.select_atoms(str_sel).indices)

# run
list_snaprange = [
    ( 0, 57877),
]

list_atompair = [
    ('h_w', 'o_1', 1),
    ('h_w', 'o_0_2', 1),
    ('h_0_1','o_w', 1),
    ('h_w','o_w', len(mda_u.select_atoms(dict_sel['o_w']))/2)
]


for tup_snaprange in list_snaprange:
    for tup_atompair in list_atompair:
        analysis.hbonds(
            universe = mda_u,
            hydrogens_sel = dict_sel[tup_atompair[0]],
            acceptors_sel = dict_sel[tup_atompair[1]],
            tup_snaprange = tup_snaprange,
            str_save = f'hbonds.{tup_atompair[0]}.{tup_atompair[1]}.{tup_snaprange[0]:07d}_{tup_snaprange[1]:07d}.csv',
            d_a_cutoff = 3.5,
            d_h_a_angle_cutoff = 150,
            int_ave = tup_atompair[2]
        )
