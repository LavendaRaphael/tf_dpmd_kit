import numpy as np

def run(
    str_data: str,
    str_save: str
) -> None:

    np_data = np.loadtxt(str_data, ndmin=2)
    np_data[:,0] += 273.15
    np.savetxt(str_save, np_data)

run('pka.csv', 'pka_kelvin.csv')
