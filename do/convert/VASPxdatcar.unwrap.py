import ovito.io
import ovito.modifiers 

ovito_pipeline = ovito.io.import_file(
    'XDATCAR',
    input_format = 'vasp',
)

ovito_pipeline.modifiers.append(ovito.modifiers.UnwrapTrajectoriesModifier())

ovito.io.export_file( 
    data = ovito_pipeline, 
    file = 'traj.lammpstrj', 
    format = 'lammps/dump',
    columns = ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"],
    multiple_frames = True
)
