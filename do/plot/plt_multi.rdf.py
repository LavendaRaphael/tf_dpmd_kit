from matplotlib import pyplot as plt
from tf_dpmd_kit import plot

str_dir_tmp = '/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_TT_H2O_126/rdf/'
dict2d_data = {
    'DPMD': {
        'h_0_1.o_w': 'rdf.h_0_1.o_w.ave.csv',
        'o_0_2.h_w': 'rdf.o_0_2.h_w.ave.csv',
        'o_1.h_w': 'rdf.o_1.h_w.ave.csv',
        'o_w.o_w': 'rdf.o_w.o_w.ave.csv',
    },
    'AIMD': {
        'h_0_1.o_w': str_dir_tmp+'rdf.h_0_1.o_w.ave.csv',
        'o_0_2.h_w': str_dir_tmp+'rdf.o_0_2.h_w.ave.csv',
        'o_1.h_w': str_dir_tmp+'rdf.o_1.h_w.ave.csv',
        'o_w.o_w': str_dir_tmp+'rdf.o_w.o_w.ave.csv',
    },
}

dict_ylim = {
    'h_0_1.o_w': (0,3),
    'o_0_2.h_w': (0,2),
    'o_1.h_w': (0,2),
    'o_w.o_w': (0,4),
}

cm=1/2.54
plot.plt_subplots(
    dict_title = {
        'h_0_1.o_w': r'H$_T$-O$_W$',
        'o_0_2.h_w': r'O$_T$-H$_W$',
        'o_1.h_w': r'$^=$O-H$_W$',
    },
    dict2d_data = dict2d_data,
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    dict_ylim = dict_ylim,
    str_save = f'rdf.dpmd.aimd.0.svg',
    tup_size = (4.3*cm, 6*cm),
    dict_legend = {
        'o_0_2.h_w': True
    }
)
plot.plt_subplots(
    dict_title = {
        'o_w.o_w': r'O$_W$-O$_W$'
    },
    dict2d_data = dict2d_data,
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    dict_ylim = dict_ylim,
    str_save = f'rdf.dpmd.aimd.1.svg',
    tup_size = (4.3*cm, 2.5*cm),
    dict_legend = {},
)

plt.show()
