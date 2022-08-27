from ovito.io import export_file
from ovito.io import import_file

ovito_pipeline = import_file(
    'XDATCAR',
    input_format = 'vasp',
)

export_file( 
    data = ovito_pipeline, 
    file = 'traj.dump', 
    format = 'lammps/dump',
    columns = ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"],
    multiple_frames = True
)
