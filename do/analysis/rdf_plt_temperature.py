from matplotlib import pyplot as plt
from tf_dpmd_kit import analysis,plm

def gen_data(
    str_pair: str,
    str_range: str,
    dict_dir: dict,
) -> dict:
    
    dict_data = {}
    for str_label, str_dir in dict_dir.items():
        dict_data[str_label] = f'{str_dir}/rdf.{str_pair}.{str_range}.csv'

    return dict_data

def plot(
    dict_title: dict,
    dict2d_data: dict,
    dict_xlim: dict,
    dict_ylim: dict,
    dict_save: dict,
) -> None:

    for str_key, str_titile in dict_title.items():
        analysis.rdf_plt_compare(
            dict_data = dict2d_data[str_key],
            tup_xlim = dict_xlim[str_key],
            tup_ylim = dict_ylim[str_key],
            str_save = dict_save[str_key],
            str_title = dict_title[str_key],
            str_xlabel = 'r (Å)',
            str_ylabel = 'g(r)',
            dict_temperature = {
                '280K': 280,
                '290K': 290,
                '300K': 300,
                '310K': 310,
                '320K': 320,
            },
            tup_colormap = (280,320)
        )

def gen_dict2d_data(
    list_pair: list,
    str_range: str,
    dict_dir: dict,
) -> dict:

    dict2d_data = {}

    for str_pair in list_pair:
        dict2d_data[str_pair] = gen_data(
            str_pair = str_pair,
            str_range = str_range,
            dict_dir = dict_dir,
        )

    return dict2d_data

def gen_dict_save(
    list_pair: list,
) -> dict:

    dict_save = {}

    for str_pair in list_pair:
        dict_save[str_pair] = f'rdf.{str_pair}.pdf'

    return dict_save

# setup

dict_title = {
    'o_0_2.h_w': r'O$_T$-H$_W$',
    'o_1.h_w': r'$^=$O-H$_W$',
    'h_0_1.o_w': r'H$_T$-O$_W$',
    'o_w.o_w': r'O$_W$-O$_W$',
}

dict_xlim = {
    'o_0_2.h_w': (1,6),
    'o_1.h_w': (1,6),
    'h_0_1.o_w': (1,6),
    'o_w.o_w': (2,6),
}

dict_ylim = {
    'o_0_2.h_w': (0,2.5),
    'o_1.h_w': (0,2.5),
    'h_0_1.o_w': (0,3),
    'o_w.o_w': (0,4),
}

dict_dir = {
    '280K': '280K/rdf/',
    '290K': '../H2CO3_TT_H2O_126/290K/rdf/',
    '300K': '../H2CO3_TT_H2O_126/300K/rdf/',
    '310K': '310K/rdf/',
    '320K': '../H2CO3_TT_H2O_126/320K/rdf/',
}
list_pair = [
    'o_0_2.h_w',
    'o_1.h_w',
    'h_0_1.o_w',
    'o_w.o_w'
]
dict2d_data = gen_dict2d_data(
    list_pair = list_pair,
    str_range = '0100000_0400000',
    dict_dir = dict_dir
)

dict_save = gen_dict_save(
    list_pair = list_pair
)

plot(
    dict_title = dict_title,
    dict2d_data = dict2d_data,
    dict_xlim = dict_xlim,
    dict_ylim = dict_ylim,
    dict_save = dict_save,
)

plt.show() 
