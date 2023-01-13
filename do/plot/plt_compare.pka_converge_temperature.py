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

str_tmp = 'dist_vp_o_1_2_fes.pka'
fig, ax = plot.plt_compare(
    dict_data = {
        '280K': f'280K/reweight/{str_tmp}.csv',
        '290K': f'290K/reweight/{str_tmp}.csv',
        '300K': f'300K/reweight/{str_tmp}.csv',
        '310K': f'310K/reweight/{str_tmp}.csv',
        '320K': f'320K/reweight/{str_tmp}.csv',
    },
    str_xlabel = 'Time (ns)',
    str_ylabel = 'pKa',
    dict_color = dict_color,
    str_save = f'reweight.{str_tmp}.converge.pdf'
)

str_tmp = 'dist_vp_o_1_2_fes.deltag'
fig, ax = plot.plt_compare(
    dict_data = {
        '280K': f'280K/reweight/{str_tmp}.csv',
        '290K': f'290K/reweight/{str_tmp}.csv',
        '300K': f'300K/reweight/{str_tmp}.csv',
        '310K': f'310K/reweight/{str_tmp}.csv',
        '320K': f'320K/reweight/{str_tmp}.csv',
    },
    str_xlabel = 'Time (ns)',
    str_ylabel = r'$\Delta$G (kJ/mol)',
    dict_color = dict_color,
    str_save = f'reweight.{str_tmp}.converge.pdf'
)
plt.show()
