from matplotlib import pyplot as plt
from tf_dpmd_kit import analysis

str_dir_aimd = '/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_TT_H2O_126/rdf/'
str_dir_dpmd = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/04.H2CO3_TT_H2O_126/330K/rdf/'

analysis.plt_subplots(
    dict_subplot = {
        #r'H$_O$-O$_W$': (0,3),
        #r'O$_H$-H$_W$': (0,2),
        r'$^=$O-H$_W$': (0,2),
        r'O$_W$-O$_W$': (0,4),
    },
    dict_data = {
        'DPMD (400 ps)': {
            r'H$_O$-O$_W$': str_dir_dpmd+'rdf.h_0_1.o_w.0050000_0130000.csv',
            r'O$_H$-H$_W$': str_dir_dpmd+'rdf.o_0_2.h_w.0050000_0130000.csv',
            r'$^=$O-H$_W$': str_dir_dpmd+'rdf.o_1.h_w.0050000_0130000.csv',
            r'O$_W$-O$_W$': str_dir_dpmd+'rdf.o_w.o_w.0050000_0130000.csv'
        },
        'AIMD (28 ps)': {
            r'H$_O$-O$_W$': str_dir_aimd+'rdf.h_0_1.o_w.0001000_0057877.csv',
            r'O$_H$-H$_W$': str_dir_aimd+'rdf.o_0_2.h_w.0001000_0057877.csv',
            r'$^=$O-H$_W$': str_dir_aimd+'rdf.o_1.h_w.0001000_0057877.csv',
            r'O$_W$-O$_W$': str_dir_aimd+'rdf.o_w.o_w.0001000_0057877.csv',
        }
    },
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    str_save = 'rdf.aimd_1.0050000_0130000.pdf',
    #tup_size = (11,5),
    bool_legend = False
)
plt.show()
