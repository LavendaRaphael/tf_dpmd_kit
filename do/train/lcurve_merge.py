import numpy as np

def lcurve_merge(
    list_data: list,
) -> None:

    np_data = np.genfromtxt(list_data[0], names=True)
    print(np_data.dtype.names)
    float_sft = np_data[-1][0]
    for str_file in list_data[1:]:
        np_tmp = np.genfromtxt(str_file, names=True)
        np_tmp['step'] += float_sft
        np_data = np.append(np_data, np_tmp)
    np.savetxt(
        fname = 'lcurve.out',
        X = np_data,
        header = ' '.join(np_data.dtype.names)
    )

lcurve_merge(
    list_data = ['lcurve.out.0', 'lcurve.out.1']
)
