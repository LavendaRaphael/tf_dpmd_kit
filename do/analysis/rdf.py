import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import numpy as np

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

# function

def def_rdf(
    str_atomgroup_0: str,
    str_atomgroup_1: str,
    list2d_range: list[list],
    dict_atomgroup: dict,
) -> None:
    int_nbins = 200
    
    mda_rdf = InterRDF( 
        dict_atomgroup[str_atomgroup_0], 
        dict_atomgroup[str_atomgroup_1],
        nbins = int_nbins,
        range=(1.0, 6.0),
        )
    for list_range in list2d_range:
        mda_rdf.run(
            start = list_range[0],
            stop = list_range[1],
            verbose = True,
            )
        array_final = np.empty( shape=(int_nbins,2) )
        array_final[:,0] = mda_rdf.results.bins
        array_final[:,1] = mda_rdf.results.rdf
        
        np.savetxt(
            fname = f'rdf.{str_atomgroup_0}.{str_atomgroup_1}.{list_range[0]:07d}_{list_range[-1]:07d}.csv',
            X = array_final,
            delimiter = ','
        )

# run
'''
list2d_range = [
    [     0,  50000],
    [ 50000, 100000],
    [100000, 150000],
    [150000, 200000],
]
#'''
list2d_range = [
    [ 50000, 200000],
]

def_rdf(
    'c',
    'o_w',
    list2d_range,
    dict_atomgroup,
)
def_rdf(
    'o_1',
    'h_w',
    list2d_range,
    dict_atomgroup,
)
def_rdf(
    'o_0_2',
    'h_w',
    list2d_range,
    dict_atomgroup,
)
def_rdf(
    'h_0_1',
    'o_w',
    list2d_range,
    dict_atomgroup,
)
def_rdf(
    'o_w',
    'o_w',
    list2d_range,
    dict_atomgroup,
)
