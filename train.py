import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def def_plt_lcurve(
    str_save: str=None
) -> None:

    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    fig, ax = plt.subplots()

    data = np.genfromtxt("lcurve.out", names=True)
    for name in ['rmse_e_trn','rmse_f_trn']:
        ax.scatter(data['step'], data[name], label=name, s=0.5)
    ax.legend(loc='lower left')
    ax.set_xlabel('Step')
    ax.set_ylabel('Loss')
    ax.set_yscale('log')
    ax.set_xlim((None,None))
    ax.grid()
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')
