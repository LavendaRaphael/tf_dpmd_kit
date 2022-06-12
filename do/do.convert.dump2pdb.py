from tf_dpmd_kit import convert
import ase.io
import numpy

array_id = numpy.arange(0,2000000,40)
str_save = 'traj.pdb'

ase_atoms = convert.def_dump2ase(
    array_id = array_id,
    type_map = ["O", "H"])

ase.io.write(format="proteindatabank", filename=str_save, images=ase_atoms)

