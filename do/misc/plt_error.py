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

str_tmp = 'dist_vp_o_1_2_fes'
#str_tmp = 'cn_o_h_fes'
#str_tmp = 'compare'

str_tmp2 = '_pka'
run(
    dict_file = {
        str_exp+'/2019_PNAS_DanielAminov/Fig_1_kelvin.csv': '2019 PNAS',
        str_exp+'/2010_JPCA_WangXiaoguang/Sfig_3_kelvin.csv': '2010 JPCA',
        str_exp+'/2009_Science_KatrinAdamczyk/pka_kelvin.csv': '2009 Science (Deuterium)',
        #str_dir+f'/02.390K_hills/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'R(V$_P$O$_C$) 390K Beta 5',
        #str_dir+f'/02.390K_hills_beta/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'R(V$_P$O$_C$) 390K Beta 10',
        str_dir+f'/03.metad/320K.5/reweight/dist_vp_o_1_2_fes{str_tmp2}.tmp.csv': 'DPMD',
        #str_dir+f'/03.metad/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'DPMD',
        #str_dir+f'/03.metad_CT/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'DPMD 2',
        #str_dir+f'/05.dist_vp.cn_o_h/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'R(V$_P$O$_C$)+cn(O$_C$H) Beta 5',
        #str_dir+f'/05.dist_vp.cn_o_h_beta/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'R(V$_P$O$_C$)+cn(O$_C$H) Beta 10',
        #str_dir+f'/06.cn_o_h/reweight.cn_o_h_fes{str_tmp2}.csv': 'cn(O$_C$H) Beta 5',
        #str_dir+f'/06.cn_o_h_beta/reweight.cn_o_h_fes{str_tmp2}.csv': 'cn(O$_C$H) Beta 10',
    },
    str_save = f'reweight.{str_tmp}{str_tmp2}_exp.pdf',
    str_ylabel = 'pKa'
)

str_tmp2 = '_deltag'
run(
    dict_file = {
        str_exp+'/2019_PNAS_DanielAminov/Fig_1_deltag.csv': '2019 PNAS',
        str_exp+'/2010_JPCA_WangXiaoguang/Sfig_3_deltag.csv': '2010 JPCA',
        str_exp+'/2009_Science_KatrinAdamczyk/deltag.csv': '2009 Science (Deuterium)',
        #str_dir+f'/02.390K_hills/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'R(V$_P$O$_C$) 390K Beta 5',
        #str_dir+f'/02.390K_hills_beta/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'R(V$_P$O$_C$) 390K Beta 10',
        str_dir+f'/03.metad/320K.5/reweight/dist_vp_o_1_2_fes{str_tmp2}.tmp.csv': 'DPMD',
        #str_dir+f'/03.metad/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'DPMD',
        #str_dir+f'/03.metad_CT/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'DPMD 2',
        #str_dir+f'/05.dist_vp.cn_o_h/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'R(V$_P$O$_C$)+cn(O$_C$H) Beta 5',
        #str_dir+f'/05.dist_vp.cn_o_h_beta/reweight.dist_vp_o_1_2_fes{str_tmp2}.csv': 'R(V$_P$O$_C$)+cn(O$_C$H) Beta 10',
        #str_dir+f'/06.cn_o_h/reweight.cn_o_h_fes{str_tmp2}.csv': 'cn(O$_C$H) Beta 5',
        #str_dir+f'/06.cn_o_h_beta/reweight.cn_o_h_fes{str_tmp2}.csv': 'cn(O$_C$H) Beta 10',
    },
    str_save = f'reweight.{str_tmp}{str_tmp2}_exp.pdf',
    str_ylabel = r'$\Delta$G (KJ/mol)'
)

plt.show()
