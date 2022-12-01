from matplotlib import pyplot as plt
from tf_dpmd_kit import analysis

def gen_dict_data(
    str_dir: str,
    str_range: str,
) -> dict:

    return {
        r'H$_0$-O$_W$': str_dir + f'rdf.h_0.o_w.{str_range}.csv',
        r'H$_1$-O$_W$': str_dir + f'rdf.h_1.o_w.{str_range}.csv',
        r'O$_0$-H$_W$': str_dir + f'rdf.o_0.h_w.{str_range}.csv',
        r'O$_1$-H$_W$': str_dir + f'rdf.o_1.h_w.{str_range}.csv',
        r'$^=$O-H$_W$': str_dir + f'rdf.o_2.h_w.{str_range}.csv',
        r'O$_W$-O$_W$': str_dir + f'rdf.o_w.o_w.{str_range}.csv',
    }


dict_data = {
    'CT': gen_dict_data('/home/faye/research_d/202203_MDCarbonicAcid/server/04.md_wall/H2CO3_CT_H2O_126/330K/rdf/', 'CT.0050000_0200000'),
}
str_range = 'CT.0050000_0200000'

analysis.plt_subplots(
    dict_subplot = {
        r'H$_0$-O$_W$': (0,3),
        r'H$_1$-O$_W$': (0,3),
        r'$^=$O-H$_W$': (0,2),
    },
    dict_data = dict_data,
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    str_save = f'rdf.{str_range}.0.pdf',
    bool_legend = False
)
analysis.plt_subplots(
    dict_subplot = {
        r'O$_0$-H$_W$': (0,2),
        r'O$_1$-H$_W$': (0,2),
        r'O$_W$-O$_W$': (0,3.5),
    },
    dict_data = dict_data,
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    str_save = f'rdf.{str_range}.1.pdf',
    bool_legend = False
)
'''
def gen_dict_data(
    str_dir: str,
    str_range: str,
) -> dict:

    return {
        r'H$_0$-O$_W$': str_dir + f'rdf.h_0.o_w.{str_range}.csv',
        r'H$_1$-O$_W$': str_dir + f'rdf.h_1.o_w.{str_range}.csv',
        r'O$_0$-H$_W$': str_dir + f'rdf.o_0.h_w.{str_range}.csv',
        r'O$_1$-H$_W$': str_dir + f'rdf.o_1.h_w.{str_range}.csv',
        r'$^=$O-H$_W$': str_dir + f'rdf.o_2.h_w.{str_range}.csv',
        r'O$_W$-O$_W$': str_dir + f'rdf.o_w.o_w.{str_range}.csv',
    }


dict_data = {
    'CT': gen_dict_data('/home/faye/research_d/202203_MDCarbonicAcid/server/04.md_wall/H2CO3_CT_H2O_126/330K/rdf/', 'CT.0050000_0200000'),
}
str_range = 'CT.0050000_0200000'

analysis.plt_subplots(
    dict_subplot = {
        r'H$_0$-O$_W$': (0,3),
        r'H$_1$-O$_W$': (0,3),
        r'$^=$O-H$_W$': (0,2),
    },
    dict_data = dict_data,
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    str_save = f'rdf.{str_range}.0.pdf',
    bool_legend = False
)
analysis.plt_subplots(
    dict_subplot = {
        r'O$_0$-H$_W$': (0,2),
        r'O$_1$-H$_W$': (0,2),
        r'O$_W$-O$_W$': (0,3.5),
    },
    dict_data = dict_data,
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    str_save = f'rdf.{str_range}.1.pdf',
    bool_legend = False
)
#'''
'''
def gen_dict_data(
    str_dir: str,
    str_range: str,
) -> dict:

    return {
        r'H$_O$-O$_W$': str_dir + f'rdf.h_0_1.o_w.{str_range}.csv',
        r'O$_H$-H$_W$': str_dir + f'rdf.o_0_2.h_w.{str_range}.csv',
        r'$^=$O-H$_W$': str_dir + f'rdf.o_1.h_w.{str_range}.csv',
        r'O$_W$-O$_W$': str_dir + f'rdf.o_w.o_w.{str_range}.csv',
    }

dict_data = {
    #'DPMD (400 ps)': gen_dict_data('/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/04.H2CO3_TT_H2O_126/330K/rdf/', 'TT.0050000_0130000'),
    'DPMD (500 ps)': gen_dict_data('/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/04.H2CO3_CC_H2O_126/330K/rdf/', 'TT.0100000_0200000'),
    'AIMD (28 ps)': gen_dict_data('/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_TT_H2O_126/rdf/', '0001000_0057877'),
}
#str_range = 'TT.0050000_0130000'
str_range = 'TT.0100000_0200000'

analysis.plt_subplots(
    dict_subplot = {
        r'H$_O$-O$_W$': (0,3),
        r'O$_H$-H$_W$': (0,2),
    },
    dict_data = dict_data,
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    str_save = f'rdf.aimd_0.{str_range}.pdf',
)

analysis.plt_subplots(
    dict_subplot = {
        r'$^=$O-H$_W$': (0,2),
        r'O$_W$-O$_W$': (0,4),
    },
    dict_data = dict_data,
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    tup_xlim = (1,6),
    str_save = f'rdf.aimd_1.{str_range}.pdf',
    bool_legend = False
)
'''
plt.show()
