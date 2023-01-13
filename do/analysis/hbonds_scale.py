import os
import numpy as np

def run(
    str_old: dict,
    float_scale: float,
    str_save: str,
) -> None:

    np_data = np.genfromtxt(str_old, names=True)
    np_data['nhbonds'] *= float_scale
    np_data['error'] *= float_scale

    np.savetxt(
        fname = str_save,
        X = np_data,
        header = ' '.join(np_data.dtype.names)
    )

run(
    str_old = 'hbonds.h_0_1.o_w.ave.temperature.csv',
    float_scale = 0.5,
    str_save = 'hbonds.h_0_1.o_w.ave.temperature.scale.csv',
)
