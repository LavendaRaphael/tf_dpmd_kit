from tf_dpmd_kit import convert
import ase.io
import numpy

array_id = numpy.arange(0,5005,5)

ase_atoms = convert.def_dump2ase(
    array_id = array_id,
    type_map = ["O", "H", "C"])

ase.io.write( filename='XDATCAR', images=ase_atoms, format='vasp-xdatcar' )
