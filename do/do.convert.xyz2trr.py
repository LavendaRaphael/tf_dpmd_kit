import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import numpy as np

float_dt = 0.00048378
float_cell = 15.679309

# common
mda_universe = mda.Universe('traj.xyz', format="XYZ", dt=float_dt)
print(mda_universe.trajectory)
mda_universe.atoms.write('traj.gro')
'''
with mda.Writer('traj.trr', n_atoms=mda_universe.atoms.n_atoms) as w:
    for ts in mda_universe.trajectory:
        ts.dimensions = [float_cell, float_cell, float_cell, 90, 90, 90]
        w.write(mda_universe.atoms)
'''
