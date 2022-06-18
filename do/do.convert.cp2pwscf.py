import numpy
import tf_dpmd_kit.convert

np_snap = numpy.arange(
    start = 50,
    stop = 3050,
    step = 50,
    )

tf_dpmd_kit.convert.def_cp2pwscf( np_snap )

