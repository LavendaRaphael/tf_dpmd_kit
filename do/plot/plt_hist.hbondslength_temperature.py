from tf_dpmd_kit import plot
import matplotlib.pyplot as plt
import math

plot.set_rcparam()
cm = 1/2.54

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

str_tmp = 'hbondslength.o_0_2.o_w.OhWaterDistance.csv'
fig, ax = plot.plt_hist(
    dict_data = {
        '280K': f'280K/hbondslength/{str_tmp}',
        '290K': f'290K/hbondslength/{str_tmp}',
        '300K': f'300K/hbondslength/{str_tmp}',
        '310K': f'310K/hbondslength/{str_tmp}',
        '320K': f'320K/hbondslength/{str_tmp}',
    },
    str_xlabel = 'Distance (Å)',
    str_ylabel = 'Probability Density',
    #tup_xlim = (0, 50),
    dict_color = dict_color,
)
plot.add_text(
    ax,
    dict_text = {
        'TT': (0.3, 0.9)
    }
)
plot.save(
    fig,
    tup_size = (4.3*cm, 3.7*cm),
    str_save = 'hbondslength.o_0_2.o_w.OhWaterDistance.hist.temperature.svg'
)

#'''
str_tmp = 'hbondslength.o_0_2.o_w.OhWaterAngle.csv'
fig, ax = plot.plt_hist(
    dict_data = {
        '280K': f'280K/hbondslength/{str_tmp}',
        '290K': f'290K/hbondslength/{str_tmp}',
        '300K': f'300K/hbondslength/{str_tmp}',
        '310K': f'310K/hbondslength/{str_tmp}',
        '320K': f'320K/hbondslength/{str_tmp}',
        #'330K': f'330K/hbondslength/{str_tmp}',
    },
    str_xlabel = r'Angle ($\degree$)',
    str_ylabel = 'Probability Density',
    tup_xlim = (0, 50),
    #tup_ylim = (None, 60),
    dict_color = dict_color,
    float_scale = 180/math.pi
)
plot.add_text(
    ax,
    dict_text = {
        'TT': (0.3, 0.9)
    }
)
plot.save(
    fig,
    tup_size = (4.3*cm, 3.7*cm),
    str_save = 'hbondslength.o_0_2.o_w.OhWaterAngle.hist.temperature.svg'
)
#'''

plt.show()
