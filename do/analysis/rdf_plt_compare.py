from matplotlib import pyplot as plt
from tf_dpmd_kit import analysis,plm

def plot(
    dict_title: dict,
    dict2d_data: dict,
    dict_xlim: dict,
    dict_ylim: dict,
    dict_save: dict,
) -> None:

    for str_key, str_title in dict_title.items():
        analysis.plt_compare_text(
            dict_data = dict2d_data[str_key],
            tup_xlim = dict_xlim[str_key],
            tup_ylim = dict_ylim[str_key],
            str_save = dict_save[str_key],
            str_title = str_title,
            str_xlabel = 'r (Å)',
            str_ylabel = 'g(r)',
        )

def gen_dict_save(
    list_pair: list,
) -> dict:

    dict_save = {}

    for str_pair in list_pair:
        dict_save[str_pair] = f'rdf.{str_pair}.compare.pdf'

    return dict_save

def gen_data(
    str_dir: dict,
    str_pair: str,
    str_range: str,
) -> dict:
    
    return f'{str_dir}/rdf.{str_pair}.{str_range}.csv'

# setup

dict_title = {
    'h_o.o_w': r'H$_O$-O$_W$',
    'o_h.h_w': r'O$_H$-H$_W$',
    'o_c.h_w': r'$^=$O-H$_W$',
}

dict_xlim = {
    'h_o.o_w': (1,6),
    'o_h.h_w': (1,6), 
    'o_c.h_w': (1,6),
}

dict_ylim = {
    'h_o.o_w': (0,3),
    'o_h.h_w': (0,2),
    'o_c.h_w': (0,2),
}

list_pair = [
    'h_o.o_w',
    'o_h.h_w',
    'o_c.h_w',
]
dict_save = gen_dict_save(
    list_pair = list_pair
)

str_dir_TT = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/H2CO3_TT_H2O_126.wall/330K/rdf/'
str_dir_CT = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/H2CO3_CT_H2O_126.wall/330K/rdf/'
str_dir_T = '.'
dict2d_data = {
    'h_o.o_w': {
        'TT': gen_data( str_dir_TT, 'h_0_1.o_w', '0100000_0400000'),
        'CT-T': gen_data( str_dir_CT, 'h_1.o_w', '0100000_0400000'),
        'CT-C': gen_data( str_dir_CT, 'h_0.o_w', '0100000_0400000'),
        r'HCO3$^-$': gen_data( str_dir_T, 'h_0.o_w', '0100000_0400000'),
    },
    'o_h.h_w': {
        'TT': gen_data( str_dir_TT, 'o_0_2.h_w', '0100000_0400000'),
        'CT-T': gen_data( str_dir_CT, 'o_1.h_w', '0100000_0400000'),
        'CT-C': gen_data( str_dir_CT, 'o_0.h_w', '0100000_0400000'),
        r'HCO3$^-$': gen_data( str_dir_T, 'o_0.h_w', '0100000_0400000'),
    },
    'o_c.h_w': {
        'TT': gen_data( str_dir_TT, 'o_1.h_w', '0100000_0400000'),
        'CT': gen_data( str_dir_CT, 'o_2.h_w', '0100000_0400000'),
        r'HCO3$^-$': gen_data( str_dir_T, 'o_1_2.h_w', '0100000_0400000'),
    }
}

plot(
    dict_title = dict_title,
    dict2d_data = dict2d_data,
    dict_xlim = dict_xlim,
    dict_ylim = dict_ylim,
    dict_save = dict_save,
)

plt.show() 
