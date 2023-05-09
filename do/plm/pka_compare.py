from tf_dpmd_kit import plot
import matplotlib.pyplot as plt
import os

str_home = os.environ['homedir']
str_exp = os.path.join(str_home, 'research/202203_MDCarbonicAcid/record/ref')

plot.plt_error(
    dict_data = {
        '2019 PNAS': str_exp+'/2019_PNAS_DanielAminov/Fig_1_kelvin.csv',
        '2016 JPCB': str_exp+'/2016_JPCB_DinePines/pka_kelvin.csv',
        '2010 JPCA': str_exp+'/2010_JPCA_WangXiaoguang/Sfig_3_kelvin.csv',
        '2009 Science': str_exp+'/2009_Science_KatrinAdamczyk/pka_kelvin.csv',
        'DPMD 1': '../03.metad/reweight.dist_vp_o_1_2_fes_pka.csv',
        'DPMD 2': 'reweight.dist_vp_o_1_2_fes.pka.csv',
    },
    str_save = f'reweight.dist_vp_o_1_2_fes.pka.compare.pdf',
    str_xlabel = r'Temperature (K)',
    str_ylabel = 'pKa',
    tup_ylim = (3.1, None)
)

plot.plt_error(
    dict_data = {
        '2019 PNAS': str_exp+'/2019_PNAS_DanielAminov/Fig_1_deltag.csv',
        '2016 JPCB': str_exp+'/2016_JPCB_DinePines/deltag.csv',
        '2010 JPCA': str_exp+'/2010_JPCA_WangXiaoguang/Sfig_3_deltag.csv',
        '2009 Science': str_exp+'/2009_Science_KatrinAdamczyk/deltag.csv',
        'DPMD 1': '../03.metad/reweight.dist_vp_o_1_2_fes_deltag.csv',
        'DPMD 2': 'reweight.dist_vp_o_1_2_fes.deltag.csv',
    },
    str_save = f'reweight.dist_vp_o_1_2_fes.deltag.compare.pdf',
    str_xlabel = r'Temperature (K)',
    str_ylabel = r'$\Delta$G (KJ/mol)',
)

plt.show()
