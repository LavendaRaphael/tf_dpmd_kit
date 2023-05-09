from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

plot.set_rcparam()
cm = 1/2.54

tup_size = (8.6*cm, 5*cm)

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

str_tmp = 'dist_vp_o_1_2_fes_pka'
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
)
plot.save(
    fig,
    tup_size = tup_size,
    str_save = f'reweight.{str_tmp}.converge.pdf'
)
plot.save(
    fig,
    tup_size = tup_size,
    str_save = f'reweight.{str_tmp}.converge.svg'
)

str_tmp = 'dist_vp_o_1_2_fes_deltag'
fig, ax = plot.plt_compare(
    dict_data = {
        '280K': f'280K/reweight/{str_tmp}.csv',
        '290K': f'290K/reweight/{str_tmp}.csv',
        '300K': f'300K/reweight/{str_tmp}.csv',
        '310K': f'310K/reweight/{str_tmp}.csv',
        '320K': f'320K/reweight/{str_tmp}.csv',
    },
    str_xlabel = 'Time (ns)',
    str_ylabel = r'$\Delta F^0$ (kJ/mol)',
    dict_color = dict_color,
    tup_ylim = (None, 25),
)
plot.save(
    fig,
    tup_size = tup_size,
    str_save = f'reweight.{str_tmp}.converge.svg'
)
plot.save(
    fig,
    tup_size = tup_size,
    str_save = f'reweight.{str_tmp}.converge.pdf'
)


plt.show()
