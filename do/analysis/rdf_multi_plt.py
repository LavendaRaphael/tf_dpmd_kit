from matplotlib import pyplot as plt
from tf_dpmd_kit import plot

str_dir_tmp = '/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_CC_H2O_126/rdf/'
dict2d_data = {
    'DPMD': {
        'cc_h_oh.o_w': 'rdf.cc_h_oh.o_w.ave.csv',
        'cc_o_oh.h_w': 'rdf.cc_o_oh.h_w.ave.csv',
        'cc_o_c.h_w': 'rdf.cc_o_c.h_w.ave.csv',
        'o_w.o_w': 'rdf.o_w.o_w.ave.csv',
    },
    'AIMD': {
        'cc_h_oh.o_w': str_dir_tmp+'rdf.h_0_1.o_w.ave.csv',
        'cc_o_oh.h_w': str_dir_tmp+'rdf.o_0_2.h_w.ave.csv',
        'cc_o_c.h_w': str_dir_tmp+'rdf.o_1.h_w.ave.csv',
        'o_w.o_w': str_dir_tmp+'rdf.o_w.o_w.ave.csv',
    },
}

dict_ylim = {
    'cc_h_oh.o_w': (0,3),
    'cc_o_oh.h_w': (0,2),
    'cc_o_c.h_w': (0,2),
    'o_w.o_w': (0,4),
}

plot.set_rcparam()
cm=1/2.54

fig, ax = plot.plt_subplots(
    dict_title = {
        'cc_h_oh.o_w': r'H$_{OH}$-O$_W$',
        'cc_o_oh.h_w': r'O$_{OH}$-H$_W$',
        'cc_o_c.h_w': r'$^=$O-H$_W$',
    },
    dict2d_data = dict2d_data,
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    dict_ylim = dict_ylim,
    dict_legend = {
        'cc_o_oh.h_w': True
    },
    float_lw = 1,
)
plot.save(
    fig,
    str_save = f'rdf.dpmd.aimd.0',
    list_type = ['pdf', 'svg'],
    tup_size = (4.3*cm, 6*cm),
)

fig,ax = plot.plt_subplots(
    dict_title = {
        'o_w.o_w': r'O$_W$-O$_W$'
    },
    dict2d_data = dict2d_data,
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    dict_ylim = dict_ylim,
    dict_legend = {},
    float_lw = 1,
)
plot.save(
    fig,
    str_save = f'rdf.dpmd.aimd.1',
    list_type = ['pdf', 'svg'],
    tup_size = (4.3*cm, 2.5*cm),
)


plt.show()
