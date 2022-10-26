import numpy as np
from tf_dpmd_kit import plm

def run(
    str_data: str,
    str_save: str
) -> None:

    np_data = np.loadtxt(str_data, ndmin=2)
    for int_i, tup_tmp in enumerate(np_data):
        np_data[int_i, 1:] = plm.pka_to_deltag(
            float_T = tup_tmp[0],
            float_pka = tup_tmp[1:]
        )
    np.savetxt(str_save, np_data)

run('pka_kelvin.csv', 'deltag.csv')
