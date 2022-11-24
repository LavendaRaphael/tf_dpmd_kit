import MDAnalysis as mda
import numpy as np
from tf_dpmd_kit import analysis

# read structrue

mda_universe = mda.Universe('traj.lammpstrj', format="LAMMPSDUMP")
mda_universe.select_atoms("type 1").types = 'O'
mda_universe.select_atoms("type 2").types = 'H'
mda_universe.select_atoms("type 3").types = 'C'
print(mda_universe.trajectory)

#'''
dict_atomgroup = {
    'o_w': mda_universe.atoms[0:126],
    'o_0_2': mda_universe.atoms[[126,128]],
    'o_1': mda_universe.atoms[[127]],
    'h_w': mda_universe.atoms[129:381],
    'h_0_1': mda_universe.atoms[[381,382]],
    'c': mda_universe.atoms[[383]]
}
#'''
'''
dict_atomgroup = {
    'o_w': mda_universe.atoms[0:126],
    'o_0_2': mda_universe.atoms[[126,127]],
    'o_1': mda_universe.atoms[[128]],
    'h_w': mda_universe.atoms[129:383]-mda_universe.atoms[[267,382]],
    'h_0_1': mda_universe.atoms[[267,382]],
    'c': mda_universe.atoms[[383]]
}
'''

for str_label, mda_atomgroup in dict_atomgroup.items():
    print(str_label, mda_atomgroup.indices)

# run
#'''
list_snaprange = [
    (     0,  50000),
    ( 50000, 100000),
    (100000, 150000),
    (150000, 200000),
    #( 50000, 200000),
]
#'''

list_atompair = [
    ('o_1', 'h_w'),
    ('o_0_2','h_w'),
    ('h_0_1','o_w'),
    ('o_w','o_w')
]

for tup_snaprange in list_snaprange:
    for tup_atompair in list_atompair:
        analysis.rdf(
            mda_atomgroup_0 = dict_atomgroup[tup_atompair[0]],
            mda_atomgroup_1 = dict_atomgroup[tup_atompair[1]],
            tup_snaprange = tup_snaprange,
            str_save = f'rdf.{tup_atompair[0]}.{tup_atompair[1]}.TT.{tup_snaprange[0]:07d}_{tup_snaprange[1]:07d}.csv',
            tup_rrange = (1.0,6.0)
        )
