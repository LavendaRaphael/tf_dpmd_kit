from matplotlib import pyplot as plt
import numpy as np
from matplotlib import rc

def grid_plt(
    str_file: str,
    str_ylabel: str = None,
    str_save: str = None,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
) -> None:

    rc('font',**{'size':20, 'family':'sans-serif','sans-serif':['Arial']})

    fig, ax = plt.subplots()
    np_data = np.loadtxt(str_file)
    ax.plot( np_data[:,0], np_data[:,1], linewidth=2)
    
    ax.set_ylabel(str_ylabel)
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.plot(1, -36, ">k", transform=ax.get_yaxis_transform(), clip_on=False)
    ax.plot(0, 47, "^k", transform=ax.get_yaxis_transform(), clip_on=False)
    if str_save:
        fig.set_size_inches(6, 5)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)

grid_plt(
    str_file = 'dist_vp_o_1_2_fes.grid',
    tup_xlim = (None, 7.5),
    str_ylabel = 'Free Energy',
    str_save = 'dist_vp_o_1_2_fes_toc.pdf',
)
plt.show()


