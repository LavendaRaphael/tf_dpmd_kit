from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

str_tmp = 'dist_vp_o_1_2_fes.10'
plm.grid_plt(
    dict_data = {
        '280K': f'280K/sum_hills/{str_tmp}.dat',
        '290K': f'290K/sum_hills/{str_tmp}.dat',
        '300K': f'300K/sum_hills/{str_tmp}.dat',
        '310K': f'310K/sum_hills/{str_tmp}.dat',
        '320K': f'320K/sum_hills/{str_tmp}.dat',
    },
    dict_temperature = {
        '280K': 280,
        '290K': 290,
        '300K': 300,
        '310K': 310,
        '320K': 320
    },
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_ylabel = 'HILLS (kJ/mol)',
    str_save = f'sum_hills.{str_tmp}.pdf',
    bool_minus = True,
    bool_maxzero = True,
    tup_xlim = (0, 14),
    #tup_ylim = (None, 60),
    tup_colormap = (280, 320),
)
plt.show()
