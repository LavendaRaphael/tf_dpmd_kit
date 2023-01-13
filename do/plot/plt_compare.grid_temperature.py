from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

dict_color = plot.dict_color_temperature(
    tup_colormap = (280, 320),
    dict_temperature = {
        '280K': 280,
        '290K': 290,
        '300K': 300,
        '310K': 310,
        '320K': 320,
    }
)

str_tmp = 'dist_vp_o_1_2_fes.ave'
fig, ax = plot.plt_compare(
    dict_data = {
        '280K': f'280K/reweight/{str_tmp}.grid',
        '290K': f'290K/reweight/{str_tmp}.grid',
        '300K': f'300K/reweight/{str_tmp}.grid',
        '310K': f'310K/reweight/{str_tmp}.grid',
        '320K': f'320K/reweight/{str_tmp}.grid',
    },
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_ylabel = 'FES (kJ/mol)',
    bool_minzero = True,
    tup_xlim = (0, 14),
    tup_ylim = (None, 60),
    dict_color = dict_color,
    str_save = f'reweight.{str_tmp}.temperature.pdf'
)

plot.inset_img(
    fig, ax,
    dict_img = {
        '../03.metad/320K/5/reweight/023719.png': (0.05, 0, 0.3, 0.4),
        '../03.metad/320K/5/reweight/023832.png': (0.2, 0.6, 0.3, 0.4),
        '../03.metad/320K/5/reweight/025097.png': (0.4, 0.15, 0.3, 0.4),
    },
    dict_arrow = {},
    str_save = f'reweight.{str_tmp}.temperature.img.pdf'
)


plt.show()
