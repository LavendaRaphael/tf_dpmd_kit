from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

str_tmp = 'dist_vp_o_1_2_fes'
plm.grid_plt(
    dict_data = {
        '280K': f'280K/reweight/{str_tmp}.grid',
        '290K': f'290K/reweight/{str_tmp}.grid',
        '300K': f'300K/reweight/{str_tmp}.grid',
        '310K': f'310K/reweight/{str_tmp}.grid',
        '320K': f'320K/reweight/{str_tmp}.grid',
    },
    dict_temperature = {
        '280K': 280,
        '290K': 290,
        '300K': 300,
        '310K': 310,
        '320K': 320
    },
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_ylabel = 'FES (kJ/mol)',
    str_save = f'reweight.{str_tmp}.pdf',
    bool_minzero = True,
    tup_xlim = (0, 14),
    tup_ylim = (None, 60),
    tup_colormap = (280, 320),
)
plt.show()
