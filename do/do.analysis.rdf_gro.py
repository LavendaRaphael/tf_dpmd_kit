import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import numpy as np

# setup
int_nbins = 200
int_step = 10
list_id = []
list_id.append( np.array([   10000,  360000]))
#list_id.append( np.array([ 400000,  800000]))
#list_id.append( np.array([ 800000, 1200000]))
#list_id.append( np.array([1200000, 1600000]))
#list_id.append( np.array([1600000, 2000000]))

float_dt = 0.00048378

# common
mda_universe = mda.Universe('traj.gro','traj.trr', dt=float_dt)
#mda_universe.select_atoms("type 1").types = 'O'
#mda_universe.select_atoms("type 2").types = 'H'
print(mda_universe.trajectory)

atomgroup_o_w = mda_universe.atoms[0:126]
print(atomgroup_o_w)
atomgroup_o_0_2 = mda_universe.atoms[[126,128]]
print(atomgroup_o_0_2)
atomgroup_o_1 = mda_universe.atoms[[127]]
print(atomgroup_o_1)
atomgroup_h_w = mda_universe.atoms[129:381]
print(atomgroup_h_w)
atomgroup_h_0_1 = mda_universe.atoms[[381,382]]
print(atomgroup_h_0_1)
atomgroup_c = mda_universe.atoms[[383]]
print(atomgroup_c)

mda_atomgroup_0 = atomgroup_c
mda_atomgroup_1 = atomgroup_o_w
str_type='c.o_w'

mda_atomgroup_0 = atomgroup_o_1
mda_atomgroup_1 = atomgroup_h_w
str_type='o_1.h_w'

mda_atomgroup_0 = atomgroup_o_0_2
mda_atomgroup_1 = atomgroup_h_w
str_type='o_0_2.h_w'

mda_atomgroup_0 = atomgroup_h_0_1
mda_atomgroup_1 = atomgroup_o_w
str_type='h_0_1.o_w'

mda_atomgroup_0 = atomgroup_o_w
mda_atomgroup_1 = atomgroup_o_w
str_type='o_w.o_w'

mda_rdf = InterRDF( 
    mda_atomgroup_0, 
    mda_atomgroup_1,
    nbins = int_nbins,
    range=(1.0, 6.0),
    )

def gen_filenpy(
        str_type,
        npint_id,
        ):
    return f'rdf.{str_type}.{npint_id[0]:07d}_{npint_id[-1]:07d}.npy'

for npint_id in list_id:
    
    npint_index = np.array(npint_id/int_step, dtype=int)
    mda_rdf.run(
        start = npint_index[0],
        stop = npint_index[1],
        verbose = True,
        )
    array_final = np.empty( shape=(2,int_nbins) )
    array_final[0] = mda_rdf.results.bins
    array_final[1] = mda_rdf.results.rdf

    str_filenpy = gen_filenpy(
        str_type = str_type,
        npint_id = npint_id,
        )
    
    np.save(
        file = str_filenpy,
        arr = array_final,
        )
    
