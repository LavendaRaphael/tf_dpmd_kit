import MDAnalysis as mda
import numpy as np
from tf_dpmd_kit import analysis

# read structrue

mda_universe = mda.Universe('traj.lammpstrj', format="LAMMPSDUMP")
mda_universe.select_atoms("type 1").types = 'O'
mda_universe.select_atoms("type 2").types = 'H'
mda_universe.select_atoms("type 3").types = 'C'
print(mda_universe.trajectory)

dict_atomgroup = {
    'o_w': mda_universe.atoms[0:126],
    'o_0_2': mda_universe.atoms[[126,128]],
    'o_1': mda_universe.atoms[[127]],
    'h_w': mda_universe.atoms[129:381],
    'h_0_1': mda_universe.atoms[[381,382]],
    'c': mda_universe.atoms[[383]]
}

# run
'''
list_snaprange = [
    (     0,  50000),
    ( 50000, 100000),
    (100000, 150000),
    (150000, 200000),
]
#'''
'''
list_snaprange = [
    ( 50000, 150000),
]
#'''
#'''
list_snaprange = [
    ( 1000, 64000),
]
#'''

list_atompair = [
    ('c', 'o_w'),
    ('o_1', 'h_w'),
    ('o_0_2','h_w'),
    ('h_0_1','o_w'),
    ('o_w','o_w')
]

for tup_atompair in list_atompair:
    analysis.rdf(
        mda_atomgroup_0 = dict_atomgroup[tup_atompair[0]],
        mda_atomgroup_1 = dict_atomgroup[tup_atompair[1]],
        list_snaprange = list_snaprange,
        str_save = f'rdf.{tup_atompair[0]}.{tup_atompair[1]}.',
        tup_rrange = (1.0,6.0)
    )
