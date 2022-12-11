import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from matplotlib.lines import Line2D
import os

def run(
    dict_file: dict,
    str_ylabel: str,
    str_save: str = None
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})
    
    fig, ax = plt.subplots()

    for str_file, str_marker in zip(dict_file, Line2D.filled_markers):
        str_label = dict_file[str_file]
        np_data = np.loadtxt(str_file, ndmin=2)
        if np_data.shape[1] > 2:
            yerr = np_data[:,2]
        else:
            yerr = None
        ax.errorbar(np_data[:,0], np_data[:,1], yerr=yerr, linestyle=':', marker=str_marker, label=str_label, capsize=3)

    ax.set_xlabel(r'Temperature (K)')
    ax.set_ylabel(str_ylabel)
    ax.legend()
    if str_save:
        fig.set_size_inches(11, 9)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)

str_home = os.environ['homedir']
str_exp = os.path.join(str_home, 'research/202203_MDCarbonicAcid/record/ref')
str_dir = os.path.join(str_home, 'research_d/202203_MDCarbonicAcid/server/05.metad')

run(
    dict_file = {
        str_exp+'/2019_PNAS_DanielAminov/Fig_1_kelvin.csv': '2019 PNAS',
        str_exp+'/2010_JPCA_WangXiaoguang/Sfig_3_kelvin.csv': '2010 JPCA',
        str_exp+'/2009_Science_KatrinAdamczyk/pka_kelvin.csv': '2009 Science (Deuterium)',
        str_dir+f'/03.metad/reweight.dist_vp_o_1_2_fes_pka.csv': 'DPMD',
    },
    str_save = f'reweight.dist_vp_o_1_2_fes_pka_exp.pdf',
    str_ylabel = 'pKa'
)

run(
    dict_file = {
        str_exp+'/2019_PNAS_DanielAminov/Fig_1_deltag.csv': '2019 PNAS',
        str_exp+'/2010_JPCA_WangXiaoguang/Sfig_3_deltag.csv': '2010 JPCA',
        str_exp+'/2009_Science_KatrinAdamczyk/deltag.csv': '2009 Science (Deuterium)',
        str_dir+f'/03.metad/reweight.dist_vp_o_1_2_fes_deltag.csv': 'DPMD',
    },
    str_save = f'reweight.dist_vp_o_1_2_fes_deltag_exp.pdf',
    str_ylabel = r'$\Delta$G (KJ/mol)'
)

plt.show()
