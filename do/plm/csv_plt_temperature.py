from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

str_tmp = 'dist_vp_o_1_2_fes_deltag'
plm.grid_plt(
    dict_data = {
        '280K': f'280K/reweight/{str_tmp}.csv',
        '290K': f'290K/reweight/{str_tmp}.csv',
        '300K': f'300K/reweight/{str_tmp}.csv',
        '310K': f'310K/reweight/{str_tmp}.csv',
        '320K': f'320K/reweight/{str_tmp}.csv',
    },
    dict_temperature = {
        '280K': 280,
        '290K': 290,
        '300K': 300,
        '310K': 310,
        '320K': 320
    },
    str_xlabel = 'Time (ns)',
    str_ylabel = r'$\Delta$G (kJ/mol)',
    str_save = f'reweight.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (16, None),
    tup_colormap = (280, 320)
)
str_tmp = 'dist_vp_o_1_2_fes_pka'
plm.grid_plt(
    dict_data = {
        '280K': f'280K/reweight/{str_tmp}.csv',
        '290K': f'290K/reweight/{str_tmp}.csv',
        '300K': f'300K/reweight/{str_tmp}.csv',
        '310K': f'310K/reweight/{str_tmp}.csv',
        '320K': f'320K/reweight/{str_tmp}.csv',
    },
    dict_temperature = {
        '280K': 280,
        '290K': 290,
        '300K': 300,
        '310K': 310,
        '320K': 320
    },
    str_xlabel = 'Time (ns)',
    str_ylabel = 'pKa',
    str_save = f'reweight.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (3.0, None),
    tup_colormap = (280, 320)
)
plt.show()
