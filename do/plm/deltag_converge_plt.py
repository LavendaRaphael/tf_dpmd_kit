from matplotlib import rc
import numpy as np
from matplotlib import pyplot as plt

def run(
    str_file: str,
    str_xheader: str,
    list_yheader: list,
    str_xlabel: str,
    list_ylabel: list,
    str_save: str = None
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    np_data = np.genfromtxt(str_file, names=True)

    fig, ax = plt.subplots()
    ax.plot( np_data[str_xheader], np_data[list_yheader[0]], linewidth=2)
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel(list_ylabel[0])

    ax1 = ax.twinx()
    ax1.plot( np_data[str_xheader], np_data[list_yheader[1]], linewidth=2)
    ax1.set_ylabel(list_ylabel[1])

    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

run(
    list_file = [
        '290K'
        'dist_vp_c_pka.csv',
    ]
    str_xheader = 'time',
    str_yheader = ['pka', 'deltag'],
    str_xlabel = 'Time (ns)',
    list_ylabel = 'pKa',
    str_save = 'dist_vp_c_pka.pdf'
)
run(
    str_file = 'dist_vp_c_pka.csv',
    str_xheader = 'time',
    str_yheader = ['pka', 'deltag'],
    str_xlabel = 'Time (ns)',
    list_ylabel = ['pKa', r'$\Delta$G (KJ/mol)'],
    str_save = 'dist_vp_c_pka.pdf'
)

plt.show()
